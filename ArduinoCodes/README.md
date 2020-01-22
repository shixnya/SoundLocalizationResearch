# Arduino codes

These codes are used to rotate a stepper motor for HRTF measurements.
This should be used together with the LabVIEW program that generates
Golay code and record sound from the mouse head.

## Hardware requirement
* Arduino + CNC kit (https://www.amazon.com/Quimat-Arduino-Contoller-Perfectly-Compatible/dp/B06XSC52SL)

## Operation
The Arduino board and the LabVIEW software communicate through
a virtual serial port. Please keep the board connected to the computer
via USB.

When the program receives a letter from the serial port, it makes a corresponding action
and print what it did in the serial port console.
Please refer to the code for available controls.
