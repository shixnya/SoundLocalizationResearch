char r = 'r';
int en = 8; // 0v
int xend = 9; // 5v x endstop
int xdir = 5; // 0v x direction 0 is clockwise, 1 is counter clockwise
int xstep = 2; // changing x step
int stepcount = 0;

void setup() {
  // put your setup code here, to run once:
  Serial.begin(9600);
  pinMode(en, OUTPUT);
  pinMode(xend, OUTPUT);
  pinMode(xdir, OUTPUT);
  pinMode(xstep, OUTPUT);
  digitalWrite(en, LOW);
  digitalWrite(xend, HIGH);
  digitalWrite(xdir, LOW);
  digitalWrite(xstep, LOW);
}

void loop() {
  // put your main code here, to run repeatedly:
  if (Serial.available() > 0) {
    r = Serial.read();
    if (r == 'r') {
      Serial.println("One Step");
      step();
      stepcount++;
    } else if (r == 'b') { // go back to origin
      Serial.println("Back to origin");
      backstep(stepcount);
      stepcount = 0;  
    } else if (r == 't') { // fine tune
      Serial.println("Ticking");
      tick();
    } else if (r == 'o') { // turn off
      Serial.println("Turn off");
      digitalWrite(en, HIGH);
    } else if (r == 'n') { // turn on
      Serial.println("Turn on");
      digitalWrite(en, LOW);
    }
  }
  delay(200);
}

void step() {
  for(int i = 0; i < 16; i++){
    tick();
  }
}

void tick() {
  delay(10);
  digitalWrite(xstep, HIGH);
  delay(10);
  digitalWrite(xstep, LOW);
}

void fasttick() {
  delay(5);
  digitalWrite(xstep, HIGH);
  delay(5);
  digitalWrite(xstep, LOW);
}


void backstep(int count) {
  digitalWrite(xdir, HIGH);
  for (int i = 0; i < 16 * count; i++) {
    fasttick();
  }
  digitalWrite(xdir, LOW);
}
