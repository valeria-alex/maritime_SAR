
# Copyright (C) 2020  Valeria Alexandra Feraru


# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
# Contact: valeria.feraru97@gmail.com

#!usr/bin/env python 


import matplotlib.pyplot as plt
import time
import math
import random
import sys
import os
from numpy import sqrt, sin, cos
import numpy as np
from matplotlib.patches import Ellipse, Polygon
from shapely.geometry import Polygon as shape_pol
from shapely.ops import cascaded_union
from descartes import PolygonPatch


def sloped(x1, y1, x2, y2): #m
	return (y2-y1)/(x2-x1)

def intercept(x1, y1, x2, y2):
	b = y2 - sloped(x1, y1, x2, y2)*x2
	return b

def checkpline(x1, y1, x2, y2, px, py):
	y = (px*sloped(x1, y1, x2, y2) + intercept(x1, y1, x2, y2))
	return y - py



def dist(x, y, z, w):
	return (sqrt( (x-z)**2+(y-w)**2 ))


def ship_path(time, lkp, speed = 10.8):
	
	return lkp + time*speed 

def elips_eq(x, lat, lon, h, k):
	result = (lon**2)*(1-((x**2- 2*x*h + h**2)/(lat**2)))
	comma = (2*k)**2 - 4*(k**2 - result)
	if comma > 0:
		y1 = ((-2*k) + sqrt(comma))/2
		y2 = ((-2*k) - sqrt(comma))/2
	else: 
		y1 = -k
		y2 = -k
	y11 = sqrt((lon**2)*(1 - (x**2)/(lat**2)))
	if math.isnan(y11):
		y11 = sqrt(-(lon**2)*(1 - (x**2)/(lat**2)))
	y21 = - y11
	return [y1, y2]

def retursoph(length, height, deliminator, positionx, positiony):
	quarter =  length/(length/deliminator)
	soph = []

	for xoxo in np.arange(positionx - length/2, positionx + length/2 + quarter, quarter): 
		eny, toy =  elips_eq(xoxo, length/2, height/2, positionx, positiony)
		soph.append([xoxo, (abs(eny))])

	for xoxo in np.arange(positionx + length/2, positionx - length/2, - quarter):
		eny, toy =  elips_eq(xoxo, length/2, height/2, positionx, positiony)
		soph.append([xoxo, (abs(toy))])

	return soph 


def normalfunction(errora, dw = False): 
	x = random.uniform(0, 100)
	if x < 71: 
		return (random.uniform(-errora, errora))
	elif x > 69 and x < 99:
		if dw == False:
			return (random.uniform(-errora*2, errora*2))
		else:
			return (random.uniform(-errora*1.5, errora*2))
	else:
		if dw == False:
			return (random.uniform(-errora*3, errora*3))
		else: 
			return (random.uniform(-errora*1.5, errora*3))


def leeway(human = False): 

	######
	wind = 8.0 #m/s -- North Sea: 1.5 to 15 m/s
	ad = 0.96/100 # a = leeway to wind ratio 0.96% wind; in some positions the ratio = 1.71 (survival suit PIW or sitting position)
	acp = 0.54/100 
	acn = -0.54/100
	#the current should be of around 0.12 - 0.25 m/s - rare cases 0.7
	bd = 0.5#7 ####THIS MAKES IT MORE ELONGATED  - total leeway speed with 0.0 we move 0.144 m (14 cm) per second; with 0.2 - 0.344 m/s = extreme situations? 
	#above can be 1 m/s
	bcp = 0.0 #CURRENTLY DOESNT WORK IF YOU CHANGE BCP OR BCN##!!!! --> should implement a function to check if a point is below or above the mid line for path planning
	bcn = 0.0
	ec = 0.094
	ed = 0.12
	#eq 

	ld = (ad)*wind + bd# +0.5*ed
	lcn = (acn)*wind + bcn 
	lcp = (acp)*wind + bcp 
	if human == True: 
		ld = (ad)*wind + bd
	return [ld, lcp + lcn] #lcp + lcn = 0

def leespeed(ld, lc, wang):
	leew_magnitude = sqrt(ld**2 + lc**2)
	new_angle = np.arcsin(lc/leew_magnitude)
	new_angle = wang - new_angle 
	return [leew_magnitude, new_angle]
	#or
def next_pos(x, y, magnitude, angle):
	x += cos(angle)*magnitude 
	y += sin(angle)*magnitude
	return [x, y] 
 
def move_elps(t, location1, location2): #THIS AND ABOVE SHOULD MOVE THE Search Area ACCORDING TO SOME ANGLE, currently not used (for simplicity)
	ld, lc = leeway()
	xled, yled = 0, 0
	print('ld, lc', ld, lc)
	magnitude, angle = leespeed(ld, lc, wang)
	
	xled, yled = next_pos(xled, yled, magnitude, angle)
	location1 = xled#0.5*t * cos(0.54) #xled * 0.01# * t*0.01
	location2 = yled
	return [location1, location2]

def piw_movement(x, y, location1, location2, t): ####this function is completly redundant now that i have the error fixed

	x = location1 #+ np.random.normal(0, 0.12)#*t  # - (1.2*t), location1 + (1.2*t))
	y = location2 #+ np.random.normal(0, 0.094)#*t  # -(0.94*t), location2 + (0.94*t))
	return [x, y]


def genpoint(a, b, k, h, angle, errors, t): ####GENERATE PIW FOR THE AREA OF THE ELLIPSE
	listpiw = []
	listaux = []
	cos_angle = np.cos(np.radians(180.-angle))
	sin_angle = np.sin(np.radians(180.-angle))
	for i in np.arange(a - k, a + k, 0.1):
		for j in np.arange(b - h, b + h, 0.1):
			xc = i - a
			yc = j - b
			xct = xc * cos_angle - yc * sin_angle
			yct = xc * sin_angle + yc * cos_angle 
			p = (xct**2.0/(k/2.0)**2.0) + (yct**2.0/(h/2.0)**2)  ### IF p < 1.0 the point should be inside the ellipse ###
			if p < 1.0:
				listaux.append([i,j])
	
	kuzco = len(listaux)
	for i in range(0, 200):
		listpiw.append([a + errors[i][0]*t, b + errors[i][1]*t])


	return listpiw
			
def insidecheck(i, j, a, b, k, h, angle, freeze, dancing): ####CHECK IF PIWs in ELLIPSE
	xc = i - a
	yc = j - b

	cos_angle = np.cos(np.radians(180.-angle))
	sin_angle = np.sin(np.radians(180.-angle))
	xct = xc * cos_angle - yc * sin_angle
	yct = xc * sin_angle + yc * cos_angle 
	p = (xct**2/(k/2.0)**2) + (yct**2/(h/2.0)**2)
	if p < 1.0:
		plt.plot(i, j,'k*')
		freeze += 1
	else:
		plt.plot(i, j, 'r*')
		dancing +=1

	return [freeze, dancing]

		

def main():
	j = 1
	fig = plt.figure()
	ax = fig.add_subplot(1, 1, 1)
	mng = plt.get_current_fig_manager()
	mng.full_screen_toggle()
	cleo = True
	listpiw = []

	dancing = 0
	freeze = 0
	downwind_leeway, crosswind_leeway = leeway()
	errord = 0.12
	errorc = 0.094
	comprised = 0
	fov = 20
	plt.ion()
	plt.grid()
	soph = [] #holds all search areas - for all timestamps
	foundline = False
	eliminate_bug = True

	elly = False

	###INITIAL POSITION###
	op = 50
	gg = 1000
	errors = []
	for error in range(200): 
		errors.append([normalfunction(0.12, False), normalfunction(0.094, False)]) ###### compute errors to simulate PIW uncertainty

	tin = 0
	ship_vel = 10.8 #m/s
	mob_signal_time = 180#seconds
	distancetolkp = ship_vel * mob_signal_time
	ax.set_facecolor('#e6f5ff')
	locationuav = [50, 1000 + distancetolkp]
	locationuav_no2 = [50, 1000 + distancetolkp]
	start_point_uav = [50, 1000 + distancetolkp]
	goals = 1
	uavspeed = 15.0 #make sure that the loops align - search for "LOOP_PART" - there should be 3: 0, 1, 2 LOOP_PART0
	sum_path = 0
	theirgoal = 1
	p = [0, 0]
	piwlkp = [50, 1000]
	locationstoplotplsnohate = []
	piw_atmob = [piwlkp[0] + downwind_leeway*(mob_signal_time) - errord*2*(mob_signal_time), piwlkp[1]]
	center_piw_atmob = [piwlkp[0] + downwind_leeway*(mob_signal_time), piwlkp[1]]

	timestamp = 10 #1 timestamp = 10 seconds

	for i in range(0, 2000, timestamp): ##########Caution: TIMESTAMPS############## LOOP_PART1
			t = i
			timothy = (time.time() - person_down)#*10 + 100

			ax = plt.gca() 
			ax.cla() # optionally clear axes

			downwind_leeway, crosswind_leeway = leeway()  

			gg = 50 + downwind_leeway*t
			op = 1000 +  crosswind_leeway*t

			exp1 = 0.12*2*2*t#**2 #############EXPAND THE AREA every second##### on X
			exp2 = 0.094*2*2*t#**2  #### ON Y ####
			theirupdown = 0
			theirlist = []
			thier_sum_path = 0
			ells = Ellipse((gg, op), exp1, exp2, el_ang, color = "#66c2ff")
			plt.plot(gg, op, color='g', marker= 'o')
			ells.set_alpha(0.4)
			ells.set_clip_box(ax.bbox)
			#ells.set_alpha(0.1)
			ax.add_artist(ells)


			if t >= mob_signal_time:
				polygon1 = []

				if foundline == False: ###FIND the intercept - where the UAV should fly to intersect with the last position where the PIW may be
					for totot in range(t, 20000):
						p = [piw_atmob[0] + downwind_leeway*(totot - mob_signal_time) - errord*2*(totot - mob_signal_time), piw_atmob[1]]

						cant = dist(locationuav[0], locationuav[1], p[0], p[1])

						auxil = p[0] - piw_atmob[0]
						if int(dist(locationuav[0], locationuav[1], p[0], p[1])/uavspeed) <= int(auxil/(downwind_leeway - 2*errord)):							
							firstone = totot - mob_signal_time
							foundline = True
							break
				pathuav = []
				pathuav.append([50, start_point_uav[1]])
				pathuav.append(p)
				plt.plot((50, p[0]),(start_point_uav[1], p[1]), 'b-')
				plt.plot(piw_atmob[0], piw_atmob[1], color ='m', marker = 'D')
				plt.plot(p[0], p[1], color ='m', marker = 'D')

				if t <= mob_signal_time + int(dist(start_point_uav[0], start_point_uav[1], p[0], p[1])/uavspeed):
					errmm = (mob_signal_time + int(dist(start_point_uav[0], start_point_uav[1], p[0], p[1])/uavspeed))
					deliminator = 0.1
				
					gg = 50 + downwind_leeway*errmm
					op = 1000 +  crosswind_leeway*errmm
		
					exp1 = ((0.12)*2)*2*errmm
					exp2 = ((0.094)*2)*2*errmm 

					  ####PLOT the search area####
					ells = Ellipse((gg, op), exp1, exp2, el_ang, color = "#ff0066")
					plt.plot(gg, op, color='g', marker= 'o')
					ells.set_alpha(0.3)
					ells.set_clip_box(ax.bbox)
					#ells.set_alpha(0.1)
					ax.add_artist(ells)

					if retursoph(exp1, exp2, deliminator, gg, op) not in soph: #merge search areas over different timestamps
					
						soph.append( retursoph(exp1, exp2, deliminator, gg, op) ) ####add polygons for each area(t)

					xs, ys = [], []
			
					for poop in range(0, len(soph[len(soph) - 1])):
						
						xs.append(round(soph[len(soph)-1][poop][0], 2)) 
						ys.append(round(soph[len(soph)-1][poop][1], 2)) 


					polygon1 = shape_pol(list(zip(xs, ys))) ### polygon in shapely (library) format ### --- so that we can merge them


					i += 1

				else:
					for i in range(len(soph) + mob_signal_time + int(dist(start_point_uav[0], start_point_uav[1], p[0], p[1])/uavspeed), t):
						deliminator = 0.1
					
						gg = 50 + downwind_leeway*i ###FIXED LEEWAY SPEED 
						op = 1000 +  crosswind_leeway*i
			
						exp1 = ((0.12)*2)*2*i#**2 
						exp2 = ((0.094)*2)*2*i ###############EXPAND THE search area #############

						if i == t:  ####PLOT the latest search area####
							ells = Ellipse((gg, op), exp1, exp2, el_ang, color = "#ff0066")
							plt.plot(gg, op, color='g', marker= 'o')
							ells.set_alpha(0.3)
							ells.set_clip_box(ax.bbox)
							#ells.set_alpha(0.1)
							ax.add_artist(ells)


						if i > 0:
							if retursoph(exp1, exp2, deliminator, gg, op) not in soph:
								soph.append( retursoph(exp1, exp2, deliminator, gg, op) ) ####add polygons for each area(t)

							xs, ys = [], []
					
							for poop in range(0, len(soph[len(soph) - 1])):
								
								xs.append(round(soph[len(soph)-1][poop][0], 2)) 
								ys.append(round(soph[len(soph)-1][poop][1], 2)) 


							polygon1 = shape_pol(list(zip(xs, ys))) ### polygon in shapely (library) format ### --- so that we can merge them


						i += 1



						
				
				if polygon1 != []:# > 0:
					if elly == True: ####PUT POLYGONS TOGETHER
						boundarii = cascaded_union([boundarii,  polygon1])
					else: 
						boundarii = cascaded_union([polygon1])
						elly = True
					patch2b = PolygonPatch(boundarii, fc='#3399ff', ec='#ff9900', alpha=1, zorder=2)
					patch2b.set_alpha(0.3)
			
				ax.add_patch(patch2b)
				plt.grid(color = '#ffcce0', linestyle=':', linewidth=2)
				




	#####################STARTING TO PLOT THE PATH#######################################
				boing = list(boundarii.exterior.coords) ###ASSIGN ALL POINTS IN FRAME POLYGON = list containing all points belonging to the polygon margins

				murs = len(boing)
				mid = boundarii.bounds[1] + (boundarii.bounds[3] - boundarii.bounds[1])/2
				upper = []
				lower = []
				########SORT ABOUT MID 
				for y in range(0, len(boing)):
					if boing[y][1] > mid: 
						upper.append(boing[y])
						#plt.plot(boing[y][0], boing[y][1], color = 'r', marker = 'o')
					elif boing[y][1] <= mid:
						lower.append(boing[y])
						#plt.plot(boing[y][0], boing[y][1], color = 'b', marker = 'o')
						#sort this??

				lower = sorted(lower)
				fov = 20

				lolo = []
				sum_path = 0

				upup = []

	########################### KEEP POINTS AT FOV DISTANCE ####################

				lolo.append(lower[0])
				for y in range(0, len(lower)):
					keep = True
					if y == 0:
						auxx = 50#
						auxy = 1000

						upper_x = 50 + downwind_leeway*t  + (0.12)*2*t#**2

						upper_y = 1000 +  crosswind_leeway*t #+ (0.094)*t

					else:
						if dist(auxx, auxy, lower[y][0], lower[y][1]) > fov:
						
							auxx = lower[y][0]
							auxy = lower[y][1]
							lolo.append(lower[y])

							if  dist(upper_x, upper_y, lower[y][0], lower[y][1]) < fov - fov/3: 
			
								break

				upup.append(upper[0])
				for y in range(0, len(upper)):																																																										
					keep = True
					if y == 0:

						auxx = 50
						auxy = 1000
					
						upper_x = 50 + downwind_leeway*t  + (0.12)*2*t#**2

						upper_y = 1000 +  crosswind_leeway*t #+ (0.094)*t

					else:
						if dist(auxx, auxy, upper[y][0], upper[y][1]) > fov:
					
							auxx = upper[y][0]
							auxy = upper[y][1]
							if int(upper[y][0]) != 50 and int(upper[y][1]) != 50:
								upup.append(upper[y])

							if  dist(upper_x, upper_y, upper[y][0], upper[y][1]) < fov - fov/3: 
						
								break

				updown = 0
				pathuav = []
				pathuav.append([50, start_point_uav[1]])
				pathuav.append(p)

				pathuav.append([upup[0][0], upup[0][1]])
				pathuav.append([upup[0][0], upup[0][1], 1000 +  crosswind_leeway*0])
				sum_path += dist(upup[0][0], upup[0][1], 50, 1000)
				upup.append([50 + downwind_leeway*t + (0.12)*2*t, 1000 +  crosswind_leeway*t])
				lolo.append([50 + downwind_leeway*t + (0.12)*2*t, 1000 +  crosswind_leeway*t])


				if len(lolo) < len(upup):
					memes = len(lolo)-1
				else:
					memes = len(upup)-1
				for mem in range(memes):
					plt.plot((upup[mem][0], lolo[mem][0]), (upup[mem][1], lolo[mem][1] ),'r-')
					sum_path += dist(upup[mem][0], upup[mem][1], lolo[mem][0], lolo[mem][1])
					if updown == 0:
			
						pathuav.append([upup[mem][0], upup[mem][1]])
	
						pathuav.append([lolo[mem][0], lolo[mem][1]])

						pathuav.append([lolo[mem+1][0], lolo[mem+1][1]])
		
						updown = 1
						sum_path += dist(lolo[mem + 1][0], lolo[mem + 1][1], lolo[mem][0], lolo[mem][1])
					else:

						updown = 0
						pathuav.append([lolo[mem][0], lolo[mem][1]])

						pathuav.append([upup[mem][0], upup[mem][1]])
						
						pathuav.append([upup[mem+1][0], upup[mem+1][1]])

				
				ey = 0
				for point in pathuav:
			
					if ey < len(pathuav) - 1:
						plt.plot((pathuav[ey + 1][0], pathuav[ey][0]), (pathuav[ey + 1][1], pathuav[ey][1] ), color='#339966', marker='_')
						ey += 1

				if t >= mob_signal_time:
					for i in range(int(uavspeed)*timestamp + 1): ######## LOOP_PART2 
						if goals < len(pathuav):
	

							p1 = locationuav
							p2 = pathuav[goals]

							x = p2[0] - p1[0]
							y = p2[1] - p1[1]
							#print (x, y)
							angle = math.atan2(y, x)

							locationuav = [p1[0] + cos(angle), p1[1] + sin(angle)]



							if int(pathuav[goals][0]) in np.arange(int(locationuav[0]) - 2.0,int(locationuav[0]) + 2.0) and int(pathuav[goals][1]) in np.arange(int(locationuav[1]) - 2.0, int(locationuav[1]) + 2.0):
								goals += 1
							plt.plot(locationuav[0], locationuav[1], color = 'k', marker = 'o') #339966

							if cleo == False:

								for point in listpiw: 
									if point[0] < locationuav[0] + fov/2 and point[0] > locationuav[0] - fov/2:
										if point[1] < locationuav[1] + fov/2 and point[1] > locationuav[1] - fov/2:
											print ('seen')
											listpiw.remove(point)



# ################ PIW ################

			if cleo and t == 0: ###########GENERATE POINTS, only once
				ggp = 50 + downwind_leeway*1 ###FIXED LEEWAY SPEED 
				opp = 1000 +  crosswind_leeway*1
	
				exp1p = ((0.12)*2)*t#**2 
				exp2p = ((0.094)*2)*t
				listpiw = genpoint(gg, op, exp1p, exp2p, el_ang, errors, t)
				cleo = False
				print ('points generated')
				for point in listpiw:
					plt.plot(point[0], point[1], color = 'r', marker = 'o')

###################MOVE POINTS #########################

			elif t > 0: #and tin < t:
				dancing = 0
				freeze = 0
				comprised = 0
		

				for point in range(len(listpiw)):
					if point/2.0 == 0:
						insteadofcrosswind_leeway = 0.1
					else:
						insteadofcrosswind_leeway = -0.1



					add = piw_movement(listpiw[point][0], listpiw[point][1], downwind_leeway*(t - tin), crosswind_leeway*(t - tin), t) #HOLD INITIAL POSITION OF POINT  
			
					listpiw[point][0] += add[0] + errors[point][0]*(t - tin) ## set list of errors created before the while loop ##
					listpiw[point][1] += add[1] + errors[point][1]*(t - tin)

		


					freeze, dancing = insidecheck(listpiw[point][0], listpiw[point][1], gg, op, exp1, exp2, el_ang, freeze, dancing) #Display PIWs: red = not within the search area
	
				tin = t
		



			prev_t = t
			if goals > 1 or t < mob_signal_time:
				plt.axis([0, 1000, 700, 1300])
			else: 
				plt.axis([0, 1000, 700, 1100 + ship_vel*mob_signal_time])


				
			

			ax.set_xlabel((t, str('s'), 'PIW spawned', 200, 'PIW remaining', len(listpiw), 'UAV takeoff (s)', mob_signal_time, 'PIW speed (cm/s)', round(downwind_leeway,2)), fontsize = 14)

			plt.draw()
			plt.pause(0.1)

	
			plt.show() 


	time.sleep(1)

#### OLD GLOBAL VARS FOR move ellipse function - not used here####
el_ang = 0.
person_down = time.time() 
wang =  0.54 #downwind is 30 deg to my position	


if __name__ == '__main__':
	try:

		main()

	except KeyboardInterrupt:
		print('Interrupted')
		try:
			sys.exit(0)
		except SystemExit:
			os._exit(0)


