# Author: Mkhanyisi Gamedze
# Benjamin Webb
# Quadruped Jumping Power Amplification
# Asymmetric Spring Loaded Inverted Pendulum model
#  April 29, 2019
# CS442 - Computational Physiology
# Team : Jacob Tower, Benjamin Webb, Mkhanyisi Gamedze

 
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import Axes3D
import math

###				 ** Global Variables **

# Initialize Global Control Variables & Parameters

# Constants
state = 0	#0 means flight, 1 means stance

dt = 0.01;	#Step size in seconds
m = 80 #mass
g = 10 #gravitaional constant
k = 500 #spring constant
restLength = 1 #Spring Rest Length

#Parameters
xvel = 5;
yvel = 0;
xpos = 0;
ypos = 1;
springXpos = 0
springYpos = 0
pivotX = 0;
xacc = 0;
yacc = 0;

dot = 0

#Variables
steps = 0;
time = 0
EPE = 0
TE = 0


# Plots the points
def updateScreenAndNextFrame(ax,yPosPlot,legAnglePlot,velocityPlot,energyPlot, dot,xpos, ypos,springXpos,springYpos,time, xvel,yvel, KE, PE, EPE):
	#Initialize
	
	# Label the figure and axes (including units)
	ax.set_xlabel("time (s)")
	ax.set_ylabel("Y0 (m)")

	yPosPlot.set_xlabel("t/s")
	yPosPlot.set_ylabel("Y0 (m)")

	legAnglePlot.set_xlabel("t/s")
	legAnglePlot.set_ylabel("Angle/Degrees")

	velocityPlot.set_xlabel("t/s")
	velocityPlot.set_ylabel("velocity/m/s")

	energyPlot.set_xlabel("t/s")
	energyPlot.set_ylabel("Mechanical Energy/J")

	# Fix axis limits so that they're the same size regardless of
	# the configuration of the arm 
	ax.set_xlim( [-1, 20] )
	ax.set_ylim( [0, 5] )

	yPosPlot.set_xlim( [0, 20] )
	yPosPlot.set_ylim( [0, 2] ) 

	#Plot
	dot, = ax.plot(xpos,ypos, 'bo')
	line, = ax.plot([xpos,springXpos], [ypos, springYpos] )

	#Plot Y position
	yPosPlot.plot(xpos, ypos, 'bo',markersize =2)


	#Plot length Angle
	angle = math.atan2( springYpos - ypos,springXpos - xpos)
	legAnglePlot.plot(time,angle,'bo',markersize =2)

	#Plot Velocity
	velocityPlot.plot(time,xvel,'bo', markersize =1, label="Vx")
	velocityPlot.plot(time,yvel,'ro',markersize =1, label = "Vy")
	blue_patch2 = mpatches.Patch(color='blue', label='Vx')
	red_patch2 = mpatches.Patch(color='red', label='Vy')
	velocityPlot.legend(handles=[blue_patch2, red_patch2])

	#Plot Energy
	energyPlot.plot(time,KE,'bo', markersize =1, label="KE")
	energyPlot.plot(time,PE,'ro',markersize =1, label = "PE")
	energyPlot.plot(time,EPE,'go', markersize =1, label="EPE")
	energyPlot.plot(time,KE+PE+EPE,'yo',markersize =1, label = "TE")
	blue_patch = mpatches.Patch(color='blue', label='KE')
	red_patch = mpatches.Patch(color='red', label='PE')
	green_patch = mpatches.Patch(color='green', label='EPE')
	yellow_patch = mpatches.Patch(color='yellow', label='TE')
	energyPlot.legend(handles=[blue_patch, red_patch, yellow_patch, green_patch])

	plt.draw()
	plt.pause(0.01)
	dot.remove()
	line.remove()

	t = 0;
	t+= dt;


	print ("Xpos, YPos:")
	print (xpos, ypos)
	print ("Xvel, Yvel:")
	print (xvel, yvel)
	print ("Xacc, Yacc:")
	print (xacc, yacc)
	print("State: ")
	print(state)
	print("\n")
	
# ASLIP equations
" Hip doesn't coincide with COM of torso, modeled as rigid body of mass 'm' and moment of inertia 'J' (J=F*d) about COM "
" ASLIP controlled by two inputs u1 (leg) and u2 (hip) "

def Df(m,J):
	return np.diag([m,m,J])

def Gf(m):
	return np.matrix([0,m*g,0])

# qs = (l, ψ, Θ)
def Ds(l, ψ, Θ):
	return np.matrix([m,0,m*l*np.cos(ψ)],[0,m*l*l,m*l*(l-l*np.sin(ψ))],[m*l*cos(ψ),m*l*(l-l*np.sin(ψ)), J+m*l*l+2*l*np.sin(ψ)])
	
# Hopping and Running	
def main3(argv):
	
	
	
	#Check number of arguments
	if (len(sys.argv) != 3):
		print ("Main function needs exactly 2 arguments: <k> <stance angle> ")
		sys.exit()
	
	# Constants
	state = 0	#0 means flight, 1 means stance

	dt = 0.01;	#Step size in seconds
	m = 80 #mass
	g = 10 #gravitaional constant
	k = 500 #spring constant
	restLength = 1 #Spring Rest Length

	#Parameters
	xvel = 5;
	yvel = 0;
	xpos = 0;
	ypos = 1;
	springXpos = 0
	springYpos = 0
	pivotX = 0;

	dot = 0

	#Variables
	steps = 0;
	time = 0
	EPE = 0
	TE = 0
	
	#Convert arguments to theta
	k = float(sys.argv[1])
	stanceAngle = np.radians(float(sys.argv[2]))

	print("Spring constant is:")
	print(k)
	print("Stance Angle (Radians) :")
	print(stanceAngle)

	
	# Plot Initialize
	fig2 = plt.figure()
	yPosPlot = fig2.add_subplot(411)

	legAnglePlot = fig2.add_subplot(412)
	velocityPlot = fig2.add_subplot(413)
	energyPlot = fig2.add_subplot(414)

	fig = plt.figure()
	ax = fig.add_subplot(111, aspect='equal')

	state = 0	# 0 -> flight, 1 -> stance
	
	##Main Loop
	while(True):

		springRestLengthY = restLength*np.cos(stanceAngle)
		springRestLengthX = restLength*np.sin(stanceAngle)

		#Figure out if Change is State has occoured
		if (state==0):
			#See if it is falling and Spring is touching the ground
			if ( (yvel < 0) and (ypos <= springRestLengthY) ):
				pivotX = xpos + springRestLengthX
				steps = steps+ 1
				state = 1

		else:
			#See if Spring Left the ground
			massPosition = np.array((xpos ,ypos))
			pivotPosition = np.array((pivotX ,0))

			compressedLength = np.absolute(np.linalg.norm(massPosition- pivotPosition))
			# print("compressedLength")
			# print(compressedLength)
			# print("pivotX")
			# print(pivotX)
			if (compressedLength > restLength):
				state = 0



		#Flight or Stance
		if (state == 0):			#Flight

			xacc = 0
			yacc = -1 * g

			#When Falling
			if (yvel<0):
				springXpos = xpos + springRestLengthX
				springYpos = ypos - springRestLengthY
			else: #Rising
				springXpos = xpos - springRestLengthX
				springYpos = ypos - springRestLengthY

			EPE = 0

		else:		#Stance
			massPosition = np.array((xpos ,ypos))
			pivotPosition = np.array((pivotX ,0))

			compressedLength = np.absolute(np.linalg.norm(massPosition- pivotPosition))

			compressedLengthsq = (xpos - pivotX)*(xpos - pivotX) + (ypos - 0)* (ypos - 0)
			compressedLength = np.sqrt(compressedLengthsq)

			compression = restLength - compressedLength

			springAngle = math.atan2( xpos - pivotX, ypos - 0)

			springForceX = k*compression* np.sin(springAngle)
			springForceY = np.absolute(k*compression*np.cos(springAngle))

			#Direction of the Force
			if (pivotX>xpos):
				springForceX = (-1)*np.absolute(springForceX)
			else:
				springForceX = np.absolute(springForceX)


			xacc = springForceX/m
			yacc = (-1*g) + ( springForceY/m)

			springXpos = pivotX
			springYpos = 0
			EPE = 0.5*k*compression*compression


		#Update Velocity and Position
		time+=dt
		xvel = xvel + xacc*dt
		yvel = yvel + yacc*dt

		xpos = xpos + xvel*dt
		ypos = ypos + yvel*dt

		totalVel = np.sqrt(xvel*xvel + yvel*yvel)
		KE = 0.5*m*(xvel*xvel + yvel*yvel)
		PE = m*g*ypos
		TE = KE + PE + EPE

		updateScreenAndNextFrame(ax,yPosPlot,legAnglePlot,velocityPlot,energyPlot, dot,xpos,ypos,springXpos,springYpos,time,xvel,yvel,KE,PE,EPE)

		print("PositionX, PositionY")
		print(xpos,ypos)
		print("VelocityX, VelocityY")
		print(xvel,yvel)
		print("Kinetic Energy")
		print(KE)
		print("Potential Energy")
		print(PE)
		print("EP Energy")
		print(EPE)
		print("Total Energy")
		print(TE)

		if (ypos < 0):
			print("Hit The Ground!!!")
			print("Numer of Steps before Failure:")
			print(steps)
			break

# Secondary Option
# Methods Developed by Ben Webb
def main2(argv):
	''' Main control loop. First creates all environmental variables and then loops through each frame and calculates
		the state of the system given change of the system since the last frame '''
		
	"""
		Useful information
		gears: (5:14:40 teeth) -> total mechanical advantage of 8
		
	"""		
		

	# Global Constant Variables
	frame=0 # each frame that is simulated
	step = 0 # each apex reached
	timer = 0.0 # total simulated time
	dt = 0.01 # each time step
	ut = 0.5 #	 update timer for the graph (ut >= dt)
	g = -9.8 # gravitational constant
	
	#Check number of arguments
	# k - spring constant of elastic material
	# stance angle - crouch/stance angle of loaded position
	# lever ratio - ratio of in lever out lever for mechanical advantage calc (directly related to ground reaction force & spring input force)
	if (len(sys.argv) != 4):
		print ("Main function needs exactly 3 arguments: <k> <stance angle> <lever ratio>")
		sys.exit()

	k =	 float(argv[1])
	angle = np.radians(float(argv[2]))
	alt = np.cos(np.radians(float(argv[2]))) # cos Θ 
	# k, angle, alt = 1000, np.radians(float(40)), np.cos(np.radians(float(40)))

	# Lever ratio gives us mechanical advantage. Compute leg parameters of outlever
	length = 10 # length of lever forelimb
	leverRatio= float(argv[3])
	x = 5 # length of other leg component (top)
	w = 0
	a = 0
	mass = 40.0
	y = (1-leverRatio) * length * np.cos(angle) # output lever length, connected to ground
	# length, leverRatio, x, w, a, mass = 10, argv, 5, 0, 0, 40.0
	
	
	# *** further calculations **
	
	# unsure here (change)
	dx = np.sin(angle) * length
	dy = np.cos(angle) * length
	dh = 0
	
	# in-lever out-lever in X and Y dimensions 
	iLX = x + dx * leverRatio
	iLY = y + dy * leverRatio 
	oLX = x - dx * (1 - leverRatio) 
	oLY	 = y - dy * (1 - leverRatio)
	p = ((iLX, oLX), (iLY, oLY))
	
	# (whats this calculation?, where did you get these?) - ask Ben
	sX = iLX-2
	sY = iLY + 8
	l = np.sqrt((5 - sX) ** 2 + (10 - sY)**2)
	
	# ************ Visualize Results & Update frames *************
	
	# Visual : Plot initialize
	fig = plt.figure()
	st = fig.suptitle("Mechanical Advantage Simulation ", fontsize="x-large")
	
	d = fig.add_subplot(211) # add a subplot to the current figure
	d.set_xlim((-5,25))
	d.set_ylim((-5,25))
	# Label the figure and axes (including units)
	d.set_xlabel("x ")
	d.set_ylabel("y ")
	d.set_title(' In-lever vs Out-lever')
	
	s = fig.add_subplot(212) # add a subplot 
	s.set_title(' Time vs In-Lever-Length')
	s.set_xlabel(" time ")
	s.set_ylabel(" In Lever Length ")
	
	# control loop, goes on and on (** Please explain what's going on here)
	while (True):
		# purpose of this?
		if (oLY <= 0.05):
			ay = (leverRatio * -k * (l - np.sqrt((iLX - sX) ** 2 + (iLY - sY) ** 2)))
		else:
			ay = (mass * g)/mass # this is simply g
			x += 5 * dt
			sX += 5 * dt

		dh += ay * dt  # Update Y velocity (Riemann Sum, dt).
		y += dh * dt  # Update Y location (Riemann Sum, dt).
		sY += dh * dt
		timer += dt
		
		dx = np.sin(angle) * length
		dy = np.cos(angle) * length
		
		# change in Lever X and Y values
		iLX = x + dx * leverRatio
		iLY = y + dy * leverRatio
		
		oLX = x - dx * (1 - leverRatio)
		oLY = y - dy * (1 - leverRatio)

		if (angle >= 0): # angle always >0 though
			# gt = mass * g * (y+(leverRatio-.5)*10)
			# t = gt - leverRatio * np.sin(angle) * k * (l - np.sqrt((iLX - sX) ** 2 + (iLY - sY) ** 2))
			# a = t / (0.136 * 100 * length**2)
			t = leverRatio * np.sin(angle) * k * (l - np.sqrt((iLX - sX) ** 2 + (iLY - sY) ** 2))
			a = t / mass
			w += a * dt
			angle += w * dt

		if frame % int(ut / dt) == 0:
			# clear plot
			#d.clf()
			
			# fresh plot
			d.plot((iLX, sX), (iLY, sY), color = "blue")
			d.plot(p[0],p[1], color = "red")
			d.plot((iLX, oLX), (iLY, oLY), color = "grey")
		
			p = ((iLX, oLX), (iLY, oLY))
			
			# d.scatter(iLX, iLY, s = 20,color = "Red")
			#d.scatter(x,y, s=20, color = "Blue")

			s.scatter(timer, np.sqrt((iLX - sX)**2+ (iLY - sY)**2), s = 10)
			
			fig.tight_layout()
			
			# shift subplots down:
			st.set_y(0.95)
			fig.subplots_adjust(top=0.85)
			
			plt.pause(ut)
		frame += 1

# kinematic model
# Quadruped jump forces and resulting dynamics Computation
def mainB(argv):
	
	
	#Check number of arguments
	# k - spring constant of elastic material
	# stance angle - crouch/stance angle of loaded position
	# lever ratio - ratio of in lever out lever for mechanical advantage calc (on the forelimb leg)
	if (len(sys.argv) != 4):
		print ("Main function needs exactly 3 arguments: <k> <stance angle> <lever ratio>")
		sys.exit()
	
	# ******	 defining system variables	********
	
	k =	 float(argv[1])
	angle = np.radians(float(argv[2]))
	alt = np.cos(np.radians(float(argv[2]))) # cos Θ 

	# Lever ratio gives us mechanical advantage. Compute leg parameters of outlever
	length = 0.15 # length of lever forelimb
	leverRatio= float(argv[3])
	# lever arms
	li=leverRatio*length
	lo=(1-leverRatio)*length
	
	# system
	m = 0.7 # kg
	g=9.81 # toward earth
	
	# angles 
	β = np.radians(30) # foot ground touch angle
	α = np.radians(90) # spring forelimb angle
	
	# other variables
	φ = np.radians(30) # dl angle equation on spring angle (look at system further then modify)
	c = 0.07 # ankle joint to top of spring distance
	
	# ******** Computations ********
	
	# spring compression
	
	# distances at angle endpoints during jump, due to spring recoil
	d1= c*c + li*li + -2*li*c*np.cos(φ)
	d2= c*c + li*li + -2*li*c*np.cos(φ+np.radians(90-angle)) # 
	
	x=np.absolute(d1-d2) # spring displacement
	
	JumpHeight = 2*k*(x**2)*(leverRatio**2)*((np.cos(β))**2)/(m*g)
	
	print (" Lever ratio :", leverRatio, "	Jump Height : ", JumpHeight )
	
# Simulation: kinematic model
# Quadruped jump forces and resulting dynamics Computation 
def mainC(argv):
	
	ut = 0.0005 # pause timestep for frame
	
	#Check number of arguments
	# k - spring constant of elastic material
	# stance angle - crouch/stance angle of loaded position
	if (len(sys.argv) != 3):
		print ("Main function needs exactly 3 arguments: <k> <stance angle>")
		sys.exit()
	
	# ******	 defining system variables	********
	
	k =	 float(argv[1])
	angle = np.radians(float(argv[2]))
	alt = np.cos(np.radians(float(argv[2]))) # cos Θ 
	
	# system
	m = 0.7 # kg
	g=9.81 # toward earth
	length = 0.10 # length of lever forelimb
	
	# angles 
	β = np.radians(30) # foot ground touch angle
	α = np.radians(90) # spring forelimb angle
	
	# other variables
	φ = np.radians(30) # dl angle equation on spring angle (look at system further then modify)
	c = 0.07 # ankle joint to top of spring distance
	
	# ******** Computations	 & Visualization ********
	
	# Visualization : Plot initialize
	fig = plt.figure()
	st = fig.suptitle("Biomechanics Simulation ", fontsize=17)
	
	s = fig.add_subplot(111) # add a subplot 
	s.set_title(' Mechanical Advantage (in-lever out-lever) vs Jump Height')
	s.set_xlabel(" Mechanical Advantage (li/lo) ")
	s.set_ylabel(" Jump Height (m) ")
	s.set_xlim((0,1))
	#s.set_ylim((0,0.4))
	Xvals = []
	Yvals =[]
	for lr in range (1,700,5):
		
		# Lever ratio gives us mechanical advantage. Compute leg parameters of outlever
		
		leverRatio= float(lr/1000)
		# lever arms
		li=leverRatio*length
		lo=(1-leverRatio)*length
	
		# spring compression
	
		# distances at angle endpoints during jump, due to spring recoil
		d1= c*c + li*li + -2*li*c*np.cos(φ)
		d2= c*c + li*li + -2*li*c*np.cos(φ+np.radians(90-angle)) # 
	
		#x=np.absolute(d1-d2) # spring displacement
		x=0.05
	
		JumpHeight = 2*k*(x**2)*(leverRatio**2)*((np.cos(β))**2)/(m*g)
		
		Xvals.append(lr)
		Yvals.append(JumpHeight)
		s.set_ylim((min(Yvals),max(Yvals)))
		s.scatter(leverRatio, JumpHeight, s = 10)
			
		
		plt.pause(ut)
		
		print (" Lever ratio :", leverRatio, "	Jump Height : ", JumpHeight )
	
			
	print ("\n\n** simulation done")
	
	plt.pause(10000) # pause till simulation stopped

# Simulation: kinematic model
# Quadruped jump forces and resulting dynamics Computation
def main(argv):
	
	ut = 0.0005 # pause timestep for frame
	
	#Check number of arguments
	# k - spring constant of elastic material
	# stance angle - crouch/stance angle of loaded position
	if (len(sys.argv) != 3):
		print ("Main function needs exactly 3 arguments: <k> <stance angle>")
		sys.exit()

	# ******	 defining system variables	  ********

	k =		float(argv[1])
	angle = np.radians(float(argv[2]))
	alt = np.cos(np.radians(float(argv[2]))) # cos Θ

	# system
	m = 0.7 # kg
	g=9.81 # toward earth
	length = 0.10 # length of lever forelimb

	# angles
	β = np.radians(30) # foot ground touch angle
	α = np.radians(90) # spring forelimb angle

	# other variables
	φ = np.radians(30) # dl angle equation on spring angle (look at system further then modify)
	c = 0.07 # ankle joint to top of spring distance
	
	# ******** Computations		& Visualization ********
	
	# Visualization : Plot initialize
	fig = plt.figure()
	st = fig.suptitle("Biomechanics Simulation ", fontsize=17)
	
	s = fig.add_subplot(111) # add a subplot
	s.set_title(' Outlever Force vs Jump Height')
	s.set_xlabel(" Out-lever force (N) ")
	s.set_ylabel(" Jump Height (m) ")
	
	#s.set_ylim((0,0.4))
	Xvals = []
	Yvals =[]
	for lr in range (1,700,5):
		
		# Lever ratio gives us mechanical advantage. Compute leg parameters of outlever
		
		leverRatio= float(lr/1000)
		# lever arms
		li=leverRatio*length
		lo=(1-leverRatio)*length
		
		# spring compression
		
		# distances at angle endpoints during jump, due to spring recoil
		d1= c*c + li*li + -2*li*c*np.cos(φ)
		d2= c*c + li*li + -2*li*c*np.cos(φ+np.radians(90-angle)) #
		
		#x=np.absolute(d1-d2) # spring displacement
		x=0.05
		
		JumpHeight = 2*k*(x**2)*(leverRatio**2)*((np.cos(β))**2)/(m*g)
		
		torque = k*x*li
		
		outforce = torque/lo
		
		Xvals.append(outforce)
		Yvals.append(JumpHeight)
		s.set_ylim((min(Yvals),max(Yvals)))
		s.set_xlim((min(Xvals),max(Xvals)))
		
		s.scatter(outforce, JumpHeight, s = 10)
		
		plt.pause(ut)
		
		print (" Lever ratio :", leverRatio, "	  Jump Height : ", JumpHeight )
	
	
	print ("\n\n** simulation done")
	
	plt.pause(10000) # pause till simulation stopped

if __name__=="__main__":
	main( sys.argv )
