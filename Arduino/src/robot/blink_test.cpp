#include <Servo.h> // Include the Servo library
#include <SoftwareSerial.h>

#include <SoftwareSerial.h>
#include <Servo.h>

Servo myServo;
SoftwareSerial mySerial(2, 3); // RX, TX pins for software serial

int servoAngle;

void setup() {
  myServo.attach(9);             // Attach servo to pin 9
  mySerial.begin(9600);          // Initialize software serial communication
  mySerial.println("Enter an angle between 0 and 180:");
}

void loop() {
  if (mySerial.available() > 0) {
    String input = mySerial.readStringUntil('\n'); // Read user input
    servoAngle = input.toInt(); // Convert input to integer

    if (servoAngle >= 0 && servoAngle <= 180) {
      myServo.write(servoAngle); // Move the servo to the specified angle
      mySerial.print("Servo moved to: ");
      mySerial.println(servoAngle);
    } else {
      mySerial.println("Invalid angle! Please enter a value between 0 and 180.");
    }
  }
}
