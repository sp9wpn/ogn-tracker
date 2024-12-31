# ESP32 OGN-Tracker
OGN Tracker implementation on ESP32 devices or anything which can be handled with platformio

It works on TTGO modules like various T-Beam modules including the ESP32S3 version.

## Compile and upload
You will need platformio which takes care for getting the right compiler and upload tools

### Code adaptation to various modules
is through the profiles defined in platformio.ini - to choose the right code use **-e** option of the **pio** command.

### Install platformio
You only need the command-line part which can be installed on Ubuntu:
```
sudo apt install cirl -y
sudo apt install curl -y
curl -fsSL -o get-platformio.py https://raw.githubusercontent.com/platformio/platformio-core-installer/master/get-platformio.py
python3 get-platformio.py
```
Try to run **pio --version** command, if not working then try **~/.platformio/penv/bin/pio --version**

Note: **DO NOT** install platformio with **sudo apt-get install platformio** if you did then remove it

### Commands related to code upload
To upload code to a **T-Beam v1.0** with **sx1276** RF chip connected on the **/dev/ttyACM0** USB-serial port
```
pio run -e ttgo-sx1276-tbeam-v10 -t upload --upload-port=/dev/ttyACM0 && minicom -D /dev/ttyACM0
```
After a succesful upload the above command starts minicom to check if the device starts properly

TTGO modules appear as **/dev/ttyACM0** or **/dev/ttyUSB0** depending on which USB-UART chip they have.

