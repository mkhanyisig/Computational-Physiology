/*
 * Author: Benjamin Webb
 * CS442 - Computational Physiology
 * Team : Jacob Tower, Benjamin Webb, Mkhanyisi Gamedze
 */
#include <stdio.h> 
#include <Servo.h> 
struct ServoControl { 
    Servo servo;
    int servoPin, analogPin, minFeedback, maxFeedback;
}; 
  
typedef struct ServoControl Struct;
Servo se1, se2, se3, se4; 
Struct s1 = {se1,12,0}; 
Struct s2 = {se2,11,1}; 
Struct s3 = {se3,10,2}; 
Struct s4 = {se4,9,3};
int minDegrees, maxDegrees;

void calibrate(Struct s, int minPos, int maxPos) {
  s.servo.write(minPos);
  minDegrees = minPos;
  delay(2000);
  s.minFeedback = analogRead(s.analogPin);
//  Serial.println(analogRead(s.analogPin));
  s.servo.write(maxPos);
  maxDegrees = maxPos;
  delay(1000);
  s.maxFeedback = analogRead(s.analogPin);
//  Serial.println(analogRead(s.analogPin));
}

void Seek(Struct s, int pos)
{
  // Start the move...
  while(abs(s.servo.read() - pos) > 2){} // wait...
}
//  s1.servo.write(a); s2.servo.write(a); s3.servo.write(a); s4.servo.write(a);

void setup(){
  Serial.begin(9600);
  delay(1000);  
//  s1.servo.attach(s1.servoPin); s3.servo.attach(s3.servoPin);
  delay(10);
  s2.servo.attach(s2.servoPin); s4.servo.attach(s4.servoPin);
  int a = 15;
//  s1.servo.write(a); s3.servo.write(180-a);
  delay(10);
  s2.servo.write(180-a); s4.servo.write(a); 
  delay(1000);
//  a = 10;
//  s1.servo.write(a); s3.servo.write(180-a);
  delay(10);
//  s2.servo.write(180-a); s4.servo.write(a);
//  delay(3000);
  s2.servo.write(180); s4.servo.write(0); 
  delay(1000);
  exit(-1);
}
