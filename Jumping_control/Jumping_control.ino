/*
 * Author: Benjamin Webb
 * CS442 - Computational Physiology
 * Team : Jacob Tower, Benjamin Webb, Mkhanyisi Gamedze
 */
 
// USART (serial bus) communication library
#include "USART.h"
#include <avr/io.h>

/* Some useful constants */
#define CW       1    // clockwise direction: to step cw, increment your index into stepMask
#define CCW     -1    // counter-clockwise direction: to step ccw, decrement your index into stepMask
#define REV    513    // our stepper has 513 steps per revolution
#define SEQ      4    // our stepper's control sequence is 4 elements long

/* If we connect the stepper's control wires {Blue, Yellow, Pink, Orange}
   to PB[3:0], then we can use the following masks on PORTB to drive the 
   stepper motor. */
const uint8_t stepClear = 0xF0;
const uint8_t stepMask[4] = { 0b00001010, 0b00000110, 0b00000101, 0b00001001 };

/* Motor position, direction, and speed parameters */
float stepTime = 0.0; // ms between steps
int stepCurrent = 0;  // the motor's current position (in steps)
int stepTarget = 0;   // the motor's target position (in steps)
int stepDir = CW;     // +1=clockwise rotation, -1=counter-clockwise rotation
int stepIdx = 0;      // the current step's index: 0-3
float stepRPM = 15.0; // revolutions per minute, rated for [0, 25] RPM at 5V
float ms = 0.0;

// controls the motors during movement
ISR( TIMER0_COMPA_vect ){
  ms+=0.1;

  // take one step once every stepTime ms
  if(ms >= stepTime){
    step();
    ms = 0.0;
  }

  // if we've reached the target, then stop (by disabling the ISR)
  if(stepCurrent == stepTarget){
    TIMSK0 &= ~0x02;
  }
}

/* Convert an angular displacement in degrees to a displacement 
   in steps. */
int deg2step( float deg ){
  int steps;
  steps =  (int) (deg*(513/360));
  return steps;
}

/* Make the stepper motor take one step. */
void step( ){
  /* Move the stepper motor by one step by overwriting PORTB with
     the next stepMask. Remember to increment or decrement stepIdx,
     depending on whether you are stepping CW or CCW. */ 
  stepIdx += stepDir;
  stepIdx = (stepIdx +4) % 4;
  PORTB = stepMask[stepIdx];
  /* Keep track of current position by incrementing or decrementing
   * stepCurrent, depending on whether you are stepping CW or CCW. */
  stepCurrent += stepDir;
  stepCurrent = (stepCurrent + 513) %513;

}

void sweep( int target, float rpm, int dir ){
  char msg[128];
  
  // Set up stepper motor control parameters: target position,
  // speed, and direction of rotation (clockwise or counter-clockwise)
  stepTarget = target;
  stepRPM = rpm; 
  if( dir==CCW ){
    stepDir = CCW;  
  } else {
    // Default to clockwise rotation
    stepDir = CW;
  }
  stepTime  = 60000.0 / (stepRPM * float(REV) );

  // Print some helpful debugging info: rotational speed and target 
  // position
  sprintf( msg, "Step time:\t%d ms\n", int(stepTime)  );
  printString( msg );
  sprintf( msg, "Target step:\t%d\n", stepTarget );
  printString( msg );

  // Trigger the beginning of the sweep, the TIMER0_COMPA ISR will
  // take care of the rest.
  step();
  ms = 0;
  TIMSK0 |= 0x02;
  while( TIMSK0 & 0x02 ){ ; }

  // Print some helpful debugging info: that we reached the target
  sprintf( msg, "Reached step:\t%d\n", stepCurrent );
  printString( msg );
}

void crouch(){
  
}

int main(){
  /* Configure Timer/Counter 0 :  0.1 ms interrupts are 
   * used to time the motor's steps, thereby controlling 
   * the speed of rotation. */
  TCCR0A = 0x02;  //CTC mode
  TCCR0B = 0x02;  //prescale by 8
  OCR0A = 200;  //compare match every .1 ms
  TIMSK0 =  0x02; //enable TC0 CompA interrupt
}
