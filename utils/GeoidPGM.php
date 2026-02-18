<?php

/**
 * GeoidPGM - Read geoid undulation (height) from a GeographicLib .pgm file.
 *
 * Supports EGM84, EGM96, and EGM2008 geoid grid files in the PGM format
 * produced by GeographicLib (e.g. egm96-5.pgm, egm2008-1.pgm).
 *
 * Download grid files from:
 *   https://geographiclib.sourceforge.io/C++/doc/geoid.html#geoidinst
 *
 * ── Usage ─────────────────────────────────────────────────────────────────────
 *   $geoid = new GeoidPGM('/path/to/egm2008-1.pgm');
 *   $h = $geoid->height($lat, $lon);          // cubic (default, more accurate)
 *   $h = $geoid->height($lat, $lon, false);   // bilinear (slightly faster)
 *   $geoid->close();                          // optional: release file handle
 *
 * ── Formula ───────────────────────────────────────────────────────────────────
 *   height_metres = offset + scale × raw_pixel_value
 *
 * ── Interpolation ─────────────────────────────────────────────────────────────
 *   Bilinear : 4-point grid-cell corners
 *   Cubic    : 12-point Karney stencil (2-D Lagrange over nodes {-1, 0, 1, 2})
 *
 * Based on Charles Karney's GeographicLib (C++ class Geoid) and PyGeodesy.
 */
class GeoidPGM
{
    // ── header ───────────────────────────────────────────────────────────────
    private float $offset;     // metres offset  (from # Offset PGM comment)
    private float $scale;      // metres / pixel (from # Scale  PGM comment)
    private int   $width;      // number of columns (longitude grid points)
    private int   $gridRows;   // number of rows    (latitude  grid points)
    private float $latRes;     // degrees per row
    private float $lonRes;     // degrees per column

    // ── file handle ──────────────────────────────────────────────────────────
    /** @var resource */
    private $fh;
    private int $dataOffset;   // byte position where pixel data begins

    // ── LRU row cache ────────────────────────────────────────────────────────
    private int   $cacheSize;
    /** @var array<int, string>  rowIndex → raw binary string (width × 2 bytes) */
    private array $rowCache  = [];
    /** @var list<int> insertion-order list used for LRU eviction */
    private array $cacheKeys = [];

    // ─────────────────────────────────────────────────────────────────────────

    /**
     * @param string $pgmFile   Path to a GeographicLib .pgm geoid file.
     * @param int    $cacheSize Number of rows to keep in RAM (LRU).
     *                          16 rows covers a full cubic stencil with margin.
     *                          Raise to e.g. 1000 for batch processing.
     * @throws RuntimeException on I/O or parse errors.
     */
    public function __construct(string $pgmFile, int $cacheSize = 16)
    {
        $this->cacheSize = max(4, $cacheSize);
        $this->openAndParse($pgmFile);
    }

    public function __destruct()
    {
        $this->close();
    }

    /** Explicitly release the file handle. Safe to call more than once. */
    public function close(): void
    {
        if (isset($this->fh) && is_resource($this->fh)) {
            fclose($this->fh);
            unset($this->fh);
        }
    }

    // ── public API ───────────────────────────────────────────────────────────

    /**
     * Return the geoid undulation N (metres) at the given WGS-84 position.
     *
     * @param float $lat   Geodetic latitude  [-90 … 90]   degrees north
     * @param float $lon   Geodetic longitude [-180 … 360] degrees east (auto-wrapped)
     * @param bool  $cubic true  = cubic interpolation  (default, ≤ 0.003 m error)
     *                     false = bilinear interpolation (≤ 0.14 m error)
     * @return float  Geoid height in metres
     */
    public function height(float $lat, float $lon, bool $cubic = true): float
    {
        // Normalise longitude → [0, 360)
        $lon = fmod($lon, 360.0);
        if ($lon < 0.0) {
            $lon += 360.0;
        }

        // Fractional grid indices.
        // Origin (row=0, col=0) = (90°N, 0°E); rows increase southward.
        $fy = (90.0 - $lat) / $this->latRes;
        $fx = $lon          / $this->lonRes;

        $row = (int) floor($fy);
        $col = (int) floor($fx);

        $dy = $fy - $row;   // sub-cell southward fraction [0, 1)
        $dx = $fx - $col;   // sub-cell eastward  fraction [0, 1)

        // Clamp row so we never step off the grid (handles exact poles)
        $row = max(0, min($this->gridRows - 2, $row));

        $pixelVal = $cubic
            ? $this->interpolateCubic($row, $col, $dy, $dx)
            : $this->interpolateBilinear($row, $col, $dy, $dx);

        return $this->offset + $this->scale * $pixelVal;
    }

    /** Return grid metadata (useful for debugging). */
    public function info(): array
    {
        return [
            'width'      => $this->width,
            'height'     => $this->gridRows,
            'latRes_deg' => $this->latRes,
            'lonRes_deg' => $this->lonRes,
            'offset_m'   => $this->offset,
            'scale_m'    => $this->scale,
        ];
    }

    // ── PGM header parser ─────────────────────────────────────────────────────

    private function openAndParse(string $path): void
    {
        $fh = fopen($path, 'rb');
        if ($fh === false) {
            throw new RuntimeException("Cannot open PGM file: $path");
        }
        $this->fh = $fh;

        // ── magic number ──
        $magic = trim(fgets($fh));
        if ($magic !== 'P5') {
            throw new RuntimeException("Not a binary PGM (P5) file. Magic='$magic'");
        }

        // ── comments + dimensions ──
        $offset = null;
        $scale  = null;
        $width  = null;
        $rows   = null;

        while (true) {
            $line = fgets($fh);
            if ($line === false) {
                throw new RuntimeException('Unexpected EOF while reading PGM header');
            }
            $line = rtrim($line);

            if (str_starts_with($line, '#')) {
                // GeographicLib embeds metadata as PGM comments
                if (preg_match('/^#\s+Offset\s+([\-\d.eE+]+)/i', $line, $m)) {
                    $offset = (float) $m[1];
                } elseif (preg_match('/^#\s+Scale\s+([\-\d.eE+]+)/i', $line, $m)) {
                    $scale = (float) $m[1];
                }
                continue;
            }

            // Dimension line(s): "width height" or just "width" then "height"
            $parts = preg_split('/\s+/', trim($line), -1, PREG_SPLIT_NO_EMPTY);
            if (count($parts) >= 2 && is_numeric($parts[0]) && is_numeric($parts[1])) {
                $width = (int) $parts[0];
                $rows  = (int) $parts[1];
                break;
            }
            if (count($parts) === 1 && is_numeric($parts[0])) {
                $width   = (int) $parts[0];
                $rowLine = fgets($fh);
                if ($rowLine === false) {
                    throw new RuntimeException('Unexpected EOF reading PGM row count');
                }
                $rows = (int) trim($rowLine);
                break;
            }
        }

        // ── max-value line (must be 65535 for 16-bit PGM) ──
        $maxVal = (int) trim(fgets($fh));
        if ($maxVal !== 65535) {
            throw new RuntimeException("Expected 16-bit PGM (maxval=65535), got $maxVal");
        }

        if ($offset === null || $scale === null) {
            throw new RuntimeException('PGM header is missing # Offset or # Scale comment');
        }
        if ($width === null || $rows === null) {
            throw new RuntimeException('Could not parse width/height from PGM header');
        }

        $this->offset   = $offset;
        $this->scale    = $scale;
        $this->width    = $width;
        $this->gridRows = $rows;
        $this->latRes   = 180.0 / ($rows - 1);
        $this->lonRes   = 360.0 / $width;

        // Remember where the binary blob starts so pixel() can fseek into it
        $this->dataOffset = ftell($fh);
    }

    // ── row fetch with LRU cache ──────────────────────────────────────────────

    /**
     * Return the raw binary string for a complete grid row (width × 2 bytes).
     * Uses an LRU cache to avoid repeated fseek/fread for the same row.
     */
    private function fetchRow(int $row): string
    {
        if (isset($this->rowCache[$row])) {
            return $this->rowCache[$row];
        }

        $bytePos = $this->dataOffset + $row * $this->width * 2;
        if (fseek($this->fh, $bytePos) !== 0) {
            throw new RuntimeException("fseek failed for row $row (byte $bytePos)");
        }

        $raw = fread($this->fh, $this->width * 2);
        if ($raw === false || strlen($raw) !== $this->width * 2) {
            throw new RuntimeException("fread failed for row $row");
        }

        // LRU eviction: remove the oldest entry when cache is full
        if (count($this->rowCache) >= $this->cacheSize) {
            $evict = array_shift($this->cacheKeys);
            unset($this->rowCache[$evict]);
        }

        $this->rowCache[$row]  = $raw;
        $this->cacheKeys[]     = $row;

        return $raw;
    }

    /**
     * Return the raw uint16 pixel at (row, col) with clamping/wrapping.
     *
     * Rows are clamped  to [0, gridRows-1]  (latitude  does not wrap).
     * Cols are wrapped  around [0, width-1] (longitude wraps at ±180°).
     */
    private function pixel(int $row, int $col): int
    {
        // Clamp latitude rows
        $row = max(0, min($this->gridRows - 1, $row));
        // Wrap longitude columns
        $col = (($col % $this->width) + $this->width) % $this->width;

        $rowData = $this->fetchRow($row);

        // Manual big-endian uint16 decode — avoids unpack() overhead
        $bytePos = $col * 2;
        return (ord($rowData[$bytePos]) << 8) | ord($rowData[$bytePos + 1]);
    }

    // ── bilinear interpolation ────────────────────────────────────────────────

    /**
     * 4-point bilinear interpolation inside a grid cell.
     *
     * Cell corners (row increases southward, col eastward):
     *   NW (row,   col  )   NE (row,   col+1)   ← dy = 0
     *   SW (row+1, col  )   SE (row+1, col+1)   ← dy = 1
     */
    private function interpolateBilinear(int $row, int $col, float $dy, float $dx): float
    {
        $nw = $this->pixel($row,     $col    );
        $ne = $this->pixel($row,     $col + 1);
        $sw = $this->pixel($row + 1, $col    );
        $se = $this->pixel($row + 1, $col + 1);

        return (1.0 - $dy) * ((1.0 - $dx) * $nw + $dx * $ne)
             +        $dy  * ((1.0 - $dx) * $sw + $dx * $se);
    }

    // ── cubic interpolation ───────────────────────────────────────────────────

    /**
     * 12-point Karney cubic interpolation — exact port of GeographicLib Geoid.cpp.
     *
     * Uses a least-squares fit of a cubic polynomial to 12 surrounding grid nodes.
     * This is NOT separable tensor-product Lagrange interpolation. The weight matrix
     * was derived by Karney via constrained least-squares and is stored as precomputed
     * integer coefficients with denominator c0 = 240.
     *
     * Stencil layout (x=dx=eastward fraction, y=dy=southward fraction):
     *
     *   node  0: pixel(row-1, col  )  x=0, y=-1
     *   node  1: pixel(row-1, col+1)  x=1, y=-1
     *   node  2: pixel(row,   col-1)  x=-1,y=0
     *   node  3: pixel(row,   col  )  x=0, y=0  <- NW corner of cell
     *   node  4: pixel(row,   col+1)  x=1, y=0  <- NE corner of cell
     *   node  5: pixel(row,   col+2)  x=2, y=0
     *   node  6: pixel(row+1, col-1)  x=-1,y=1
     *   node  7: pixel(row+1, col  )  x=0, y=1  <- SW corner of cell
     *   node  8: pixel(row+1, col+1)  x=1, y=1  <- SE corner of cell
     *   node  9: pixel(row+1, col+2)  x=2, y=1
     *   node 10: pixel(row+2, col  )  x=0, y=2
     *   node 11: pixel(row+2, col+1)  x=1, y=2
     *
     * Polynomial basis (10 terms): [1, x, y, x^2, xy, y^2, x^3, x^2y, xy^2, y^3]
     *
     * weight_i(dx,dy) = (1/c0) * sum_j( c3[i][j] * basis_j(dx,dy) )
     *
     * c3[12][10] and c0=240 taken verbatim from GeographicLib/src/Geoid.cpp
     * (Charles Karney, MIT licence).
     *
     * Reference: https://geographiclib.sourceforge.io/C++/doc/geoid.html
     */
    private function interpolateCubic(int $row, int $col, float $dy, float $dx): float
    {
        // Karney weight table, denominator = 240.
        // Source: GeographicLib/src/Geoid.cpp  (const int Geoid::c3_[])
        static $C0 = 240;
        static $C3 = [
        //    1    x    y   x2   xy   y2   x3  x2y  xy2   y3
            [  9, -18, -88,   0,  96,  90,   0,   0, -60, -20],  //  0: (row-1, col  )
            [ -9,  18,   8,   0, -96,  30,   0,   0,  60, -20],  //  1: (row-1, col+1)
            [  9, -88, -18,  90,  96,   0, -20, -60,   0,   0],  //  2: (row,   col-1)
            [186, -42, -42,-150, -96,-150,  60,  60,  60,  60],  //  3: (row,   col  )
            [ 54, 162, -78,  30, -24, -90, -60,  60, -60,  60],  //  4: (row,   col+1)
            [ -9, -32,  18,  30,  24,   0,  20, -60,   0,   0],  //  5: (row,   col+2)
            [ -9,   8,  18,  30, -96,   0, -20,  60,   0,   0],  //  6: (row+1, col-1)
            [ 54, -78, 162, -90, -24,  30,  60, -60,  60, -60],  //  7: (row+1, col  )
            [-54,  78,  78,  90, 144,  90, -60, -60, -60, -60],  //  8: (row+1, col+1)
            [  9,  -8, -18, -30, -24,   0,  20,  60,   0,   0],  //  9: (row+1, col+2)
            [ -9,  18, -32,   0,  24,  30,   0,   0, -60,  20],  // 10: (row+2, col  )
            [  9, -18,  -8,   0, -24, -30,   0,   0,  60,  20],  // 11: (row+2, col+1)
        ];

        // Read the 12 stencil pixels
        $v = [
            (float) $this->pixel($row - 1, $col    ),  //  0
            (float) $this->pixel($row - 1, $col + 1),  //  1
            (float) $this->pixel($row,     $col - 1),  //  2
            (float) $this->pixel($row,     $col    ),  //  3
            (float) $this->pixel($row,     $col + 1),  //  4
            (float) $this->pixel($row,     $col + 2),  //  5
            (float) $this->pixel($row + 1, $col - 1),  //  6
            (float) $this->pixel($row + 1, $col    ),  //  7
            (float) $this->pixel($row + 1, $col + 1),  //  8
            (float) $this->pixel($row + 1, $col + 2),  //  9
            (float) $this->pixel($row + 2, $col    ),  // 10
            (float) $this->pixel($row + 2, $col + 1),  // 11
        ];

        // Polynomial basis evaluated at (dx, dy)
        $x  = $dx;
        $y  = $dy;
        $b  = [
            1.0,           // 1
            $x,            // x
            $y,            // y
            $x * $x,       // x^2
            $x * $y,       // xy
            $y * $y,       // y^2
            $x * $x * $x,  // x^3
            $x * $x * $y,  // x^2*y
            $x * $y * $y,  // x*y^2
            $y * $y * $y,  // y^3
        ];

        // Weighted sum: (1/c0) * sum_i( v[i] * sum_j( c3[i][j]*b[j] ) )
        $result = 0.0;
        for ($i = 0; $i < 12; $i++) {
            $w = 0;
            for ($j = 0; $j < 10; $j++) {
                $w += $C3[$i][$j] * $b[$j];
            }
            $result += $w * $v[$i];
        }

        return $result / $C0;
    }
}
?>
