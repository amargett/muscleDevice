#include <Servo.h> // Include the Servo library

Servo myServo; // Create a servo object to control the servo
int servoAngle; // Variable to store the user-defined servo angle (0-180 degrees)

void setup() {
  myServo.attach(9);       // Attach the servo to pin 9
  Serial.begin(9600);      // Start serial communication at 9600 baud
  Serial.println("Enter an angle between 0 and 180:"); // Prompt the user
}

void loop() {
  if (Serial.available() > 0) { // Check if data is available from the serial monitor
    String input = Serial.readStringUntil('\n'); // Read the user input until a newline
    servoAngle = input.toInt(); // Convert the input to an integer

    // Validate the input
    if (servoAngle >= 0 && servoAngle <= 180) {
      myServo.write(servoAngle); // Move the servo to the specified angle
      Serial.print("Servo moved to: ");
      Serial.println(servoAngle);
    } else {
      Serial.println("Invalid angle! Please enter a value between 0 and 180.");
    }
  }
}
