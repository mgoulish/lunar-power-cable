#! /usr/bin/env python3

import math
import sys
from PIL import Image, ImageDraw, ImageFont, ImageOps
import matplotlib.pyplot as mpl
import numpy as np
from skimage.draw import (line, polygon, disk,
                          circle_perimeter,
                          rectangle,
                          ellipse, ellipse_perimeter,
                          bezier_curve)




#------------------------------------------------
# Draw the concentric shell cross-sections at a 
# particular time step.
#------------------------------------------------
def draw ( n_time_step ) :
  # Define colors for the temperature breaks,
  # in Kelvins
  K_320   = (1.00, 1.00, 0.00)
  K_310   = (1.00, 0.75, 0.00)
  K_300   = (1.00, 0.53, 0.00)
  K_290   = (1.00, 0.25, 0.00)
  K_280   = (0.93, 0.00, 0.00)
  K_270   = (0.84, 0.00, 0.00)
  K_260   = (0.74, 0.00, 0.00)
  K_250   = (0.60, 0.00, 0.00)
  K_240   = (0.45, 0.00, 0.00)
  AMBIENT = (0.00, 0.00, 0.00)   # 230


  # Make the image.
  image_size = 1200
  center = image_size / 2
  img = np.zeros((image_size, image_size, 3), dtype=np.double)
  print ( f"drawing image for timestep {n_time_step}" )
  pixels_per_cm = 15

  # Iterate through all the concentric shells of regolith
  # around the power cable, choose the color for their temp
  # and draw them.
  # Note that we are iterating backwards here, from largest 
  # to smallest shell, so that the smaller circles will get 
  # drawn on top of the larger ones.
  for s in range(299, -1, -1) :
    temp = shell_temps[s]
    if temp <= ambient_temperature :
      continue
    if temp   >= 320 :
      color = K_320
    elif temp >= 310 :
      color = K_310
    elif temp >= 300 :
      color = K_300
    elif temp >= 290 :
      color = K_290
    elif temp >= 280 :
      color = K_280
    elif temp >= 270 :
      color = K_270
    elif temp >= 260 :
      color = K_260
    elif temp >= 250 :
      color = K_250
    elif temp >= 240 :
      color = K_240
    else :
      color = AMBIENT
    # 's' is the shell number, which is the same as its radius,
    # in centimeters. 
    r = (s+1) * pixels_per_cm
    rr, cc = disk((center, center), r, shape=img.shape)
    img[rr, cc, :] = color

  # Draw the 25 cm line ----------------------------
  # i.e. 25 cm above the cable
  y = int(center) - 15 * 25
  rr, cc = line ( y, 0, y, image_size-1)
  img[rr, cc] = ( 1.0, 1.0, 1.0 )

  # Clear matplotlib figure from previous image.
  # Or we will get new text superimposed on old text.
  mpl.clf()
  mpl.figure(figsize=(4,4))

  # Draw the temperature key --------------------------------
  rect_y      = 1100
  rect_x      = 50
  rect_width  = 120
  rect_height = 60
  temp = 240
  for color in (K_240,K_250,K_260,K_270,K_280,K_290,K_300,K_310,K_320):
    rect_start  = (rect_y, rect_x)
    rect_extent = (rect_height, rect_width)
    rr, cc = rectangle(rect_start, extent=rect_extent, shape=img.shape)
    img[rr, cc] = color
    mpl.text ( rect_x, rect_y - 5, str(temp), color=(1,1,1) )
    rect_x += rect_width
    temp += 10

  # Let's display number of days.
  seconds_per_day = 86400
  timesteps_per_day = seconds_per_day / time_step
  days = int(n_time_step / timesteps_per_day)
  label = str(days) + " Days"
  mpl.text ( 100, 100, label, color=(1,1,1) )

  # Label the 25 cm line
  mpl.text ( 900, y, "25 cm", color=(1,1,1) )

  # Save out the image
  filename = './frames/img_' + "{:05d}".format(n_time_step) + '.png'
  mpl.imshow ( img )
  mpl.axis ( 'off' )
  mpl.savefig ( filename, bbox_inches='tight', dpi=400 )



#============================================
# Main 
#============================================


#--------------------------------------------
# Initial Conditions
#--------------------------------------------
sim_duration        = 36000 # In time-steps
time_step           = 900  # Size of each time step in seconds.
# This is a little over 1 year

regolith_density    = 1.8   # g/cm3
reg_thermal_conduct = 0.85  # Watts / ( cm * K)
reg_spec_heat       = 1.512 # J/(g * K)
aluminum_density    = 2.7   # g/cm3
cable_length        = 100   # cm
cable_diameter      = 2     # cm. This is a *big* cable.
cable_mass          = math.pi * (cable_diameter/2)**2 * cable_length * \
                      aluminum_density
ambient_temperature = 230   # 230 K ~= -46 F
heat_dissipation    = 1     # Watt, or Joule per second, in the cable
al_spec_heat        = 0.903 # Joules / (gram * K)
original_cable_temp = ambient_temperature
added_cable_joules  = 0
min_temp_delta      = 0.1 

#-------------------------------------------
# find cable surface area
#-------------------------------------------
cable_area = cable_diameter * math.pi * 100
cable_area /= 10000  # Convert to square meters.
#print ( f"cable area is {cable_area} m^2" )


#---------------------------------------------------
# We will look at 1-cm-thick shells going outward
# from the cable. The cable is 2 m deep so let's 
# look at 300 of them. (The top 100 shells will be
# truncated by the surface.)
# These arrays will track -- or at least store --
# all the quantities I care about for each shell.
#---------------------------------------------------
n_shells = 300
shell_temps  = [ambient_temperature] * n_shells
shell_energy = [0] * n_shells
shell_areas  = [0] * n_shells
shell_masses = [0] * n_shells


#-----------------------------------------------
# Let's pre-calc the masses and areas of all 
# the shells. Not shell 0, though. That's the cable.
#-----------------------------------------------
shell_masses[0] = cable_mass
shell_areas [0] = cable_area
shell_energy[0] = cable_mass * ambient_temperature * al_spec_heat
for r in range(1, 300) : # r stands for radius, in cm
  shell_area      = math.pi * 2 * r * 100
  shell_areas[r]  = shell_area / 10000 # In m2
  shell_volume    = shell_area # It's 1 cm thick.
  shell_mass      = shell_volume * regolith_density # In grams
  shell_masses[r] = shell_mass
  shell_energy[r] = shell_masses[r] * ambient_temperature * reg_spec_heat
  #print ( f"Mass of shell {r} is {shell_mass}" )


#-----------------------------------------------
# Once per step for the durtion of the sim,
# propagate heat energy from the cable outward.
#-----------------------------------------------
for step in range(0, sim_duration+1) :
  if 0 == (step % 1000) :
    print ( f"step {step}" )
  for s in range(0, 299) :
    # heat diff between this shell and next one out
    temp_diff = shell_temps[s] - shell_temps[s+1]
    if temp_diff < min_temp_delta :
      break
    # Thermal conductivity of regolith at depth, 
    # which I take to mean thermal conductivity of 
    # regolith compressed a little, is:
    #   8.5e–3 W m–1 K–1 at a depth of 1 m
    # So for a centimeter instead of a meter, that is:
    #   8.5e-1 W / ( cm * K )
    # It changes inversely with length. I think.
    # And I think that means per meter squared.
    # So first find how much energy we are losing
    # per square meter because of the temp diff.
    # The TC is Watts per K, so scale for diff
    heat_loss = reg_thermal_conduct * temp_diff
    # And this is per square meter of surface, 
    # so now scale by this shell's surface area.
    heat_loss *= shell_areas[s]
    # And this heat loss is in Watts, so -- since 
    # we are running this sim in N-second time steps --
    # we multiply by N to get Joules.
    heat_loss *= time_step 
    # So -- remove this quantity of heat energy, in Joules,
    # from this shell and give it to the next one!
    shell_energy[s]   -= heat_loss
    shell_energy[s+1] += heat_loss

    # Calc new temps for these two shells
    # This shell
    joules_per_gram = shell_energy[s] / shell_masses[s]
    kelvins = joules_per_gram / reg_spec_heat
    shell_temps[s] = kelvins 
    # The next shell
    joules_per_gram = shell_energy[s+1] / shell_masses[s+1]
    kelvins = joules_per_gram / reg_spec_heat
    shell_temps[s+1] = kelvins 
    # end of heat propagation loop --------------

  #-----------------------------------------------
  # Now that we have propagated as much energy
  # as we can outward, add new energy to the
  # cable.   Shell 0 is the cable.
  #-----------------------------------------------
  shell_energy[0] += heat_dissipation * time_step 
  joules_per_gram = shell_energy[0] / shell_masses[0]
  kelvins = joules_per_gram / al_spec_heat
  shell_temps[0] = kelvins  


  #--------------------------------------------------
  # Make a drawing every 2880 timesteps.
  # 2880 * 15 minutes is 30 days.
  #--------------------------------------------------
  if step > 0 and 0 == (step % 2880) :
    draw ( step )
  # end of sim loop ----------------------
  




