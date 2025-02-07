const int magnetPin = 9; // Pin connected to the electromagnet

void setup() {
  pinMode(magnetPin, OUTPUT); // Set the electromagnet pin as output
  digitalWrite(magnetPin, LOW); // Ensure the electromagnet is off at startup
  Serial.begin(9600);          // Start serial communication
  Serial.println("Type 'ON' to turn on the electromagnet and 'OFF' to turn it off:");
}

void loop() {
  if (Serial.available() > 0) {          // Check if there is data in the Serial Monitor
    String command = Serial.readStringUntil('\n'); // Read the input
    command.trim();                      // Remove any leading/trailing whitespace

    if (command == "ON") {               // If the user types "ON"
      digitalWrite(magnetPin, HIGH);     // Turn on the electromagnet
      Serial.println("Electromagnet is ON");
    } else if (command == "OFF") {       // If the user types "OFF"
      digitalWrite(magnetPin, LOW);      // Turn off the electromagnet
      Serial.println("Electromagnet is OFF");
    } else {                             // Invalid command
      Serial.println("Invalid command! Type 'ON' or 'OFF'.");
    }
  }
}
