import serial
import time
from astropy.time import Time
from serial.tools import list_ports

COMM = 6

try:
    serial_port = serial.Serial(port=f"COM{COMM}",
                                baudrate=115200,
                                bytesize=8,
                                timeout=2,
                                stopbits=serial.STOPBITS_ONE,
                                parity=serial.PARITY_NONE
                                )
except serial.SerialException:
    print(f"Could not open serial port on COM{COMM}")
    print("Available ports:")
    ports = serial.tools.list_ports.comports()
    for port, desc, hwid in sorted(ports):
        print(f"{port}: {desc} [{hwid}]")
    exit(1)

time.sleep(2)
serial_port.write(b"\n")
logged_in = 0
sys_time = None
dvc_time = None
old_time = None

try:
    while True:
        string = serial_port.readline()
        string = string.decode("Ascii").strip()
        if string:
            print("\t" + string)
        else:
            continue

        if logged_in == 0 and "login" in string.lower():
            serial_port.write(b"\n")  # insert login
            logged_in = 1

        if logged_in == 1 and "password" in string.lower():
            serial_port.write(b"\n")  # insert password
            logged_in = 2

        if logged_in == 2:
            time.sleep(2)
            serial_port.write(b"date\n")
            logged_in = 3

        if logged_in == 3 and "UTC" in string:
            dvc_time = string
            serial_port.write(b"date\n")
            logged_in = 4

        if logged_in == 4 and "UTC" in string:
            if string != dvc_time:
                old_time = dvc_time
                dvc_time = string
                sys_time = Time.now()
                logged_in = 5
            else:
                serial_port.write(b"date\n")

        if logged_in == 5:
            serial_port.write(b"exit\n")
            break
except KeyboardInterrupt:
    serial_port.write(b"exit\n")
finally:
    print("Closing port")
    serial_port.close()


if all([sys_time, dvc_time, old_time]):
    format_string = "%a %b %d %H:%M:%S %Z %Y"

    # Parse the time string into a datetime object
    parsed_time = Time.strptime(dvc_time.strip().replace("  ", " "), format_string)
    parsed_time.format = "iso"

    parsed_old_time = Time.strptime(old_time.strip().replace("  ", " "), format_string)
    parsed_old_time.format = "iso"

    print("PROGRAM SUCCESSFUL")
    print(f"System time was {sys_time} when the device time changed from {parsed_old_time} to {parsed_time}")
    print(f"Difference (system time is ahead of device by) is {(sys_time - parsed_time).to("s").value:.3f} seconds")
