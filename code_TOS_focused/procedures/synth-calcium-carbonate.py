##########################################################################################################
# Calcium Carbonate Polymorphs
##########################################################################################################

###########################
#    Import Libraries     #
###########################

# in-house Chemputer libraries and related
from chempiler import Chempiler
import xdl
import ChemputerAPI
import commanduinolabware
import logging

# for file handling
import appdirs
import os
import datetime

# for automatic emailing
import smtplib, ssl
import getpass

# for data manipulation
import pandas as pd



###########################
#    Chempiler Object     #
###########################

# set working directory to same as this script
os.chdir(os.path.dirname(os.path.abspath(__file__)))
path=os.getcwd()

# define graph and xdl files
GRAPH_FILE = path + "general-graph.json"

# create Chempiler object
c = Chempiler(
    experiment_code = "calcium-carbonate-polymorphs", 
    output_dir=path, 
    simulation=True, 
    graph_file=GRAPH_FILE, 
    device_modules=[ChemputerAPI, commanduinolabware]
)



###########################
#   User-Set Variables    #
###########################

limiting_reagent_moles = 0.02 #0.04   # mol, limiting reagent: CaCl2.2H2O
average_wash = 0   # mL



###########################
#   Automatic Variables   #
###########################

#define basic calculations for amounts of substance
def calculate_mass(reagent, equiv) -> float:
    MW = reagent_MW.get(str(reagent))
    mass = float(MW) * (float(limiting_reagent_moles) * float(equiv))
    return round(mass, 2)

def calculate_volume_soln(reagent, equiv) -> float:
    conc = reagent_conc.get(str(reagent))
    volume_L = (float(limiting_reagent_moles) * float(equiv)) / float(conc)
    volume_mL = volume_L * 1000
    return round(volume_mL, 2)

def calculate_volume_pure(reagent, conc_wrt_LR) -> float:   
    volume_L = limiting_reagent_moles / conc_wrt_LR   # L
    volume_mL = volume_L * 1000
    return round(volume_mL, 2)

#define dictionaries containing physical data on reagents
reagent_MW = {
    "CaCl2.2H2O":"147.01",
    "K2CO3":"138.21",
    "water":"18.02", 
}   # g/mol
reagent_conc = {
    "CaCl2.2H2O":"2",
    "K2CO3":"2"
}   # Molar
reagent_density = {
    "water":"1.000"
}   # g/mL

#calculate amounts of substances
amount_CaCl2 = calculate_volume_soln("CaCl2.2H2O", 1)   # mL
amount_K2CO3 = calculate_volume_soln("K2CO3", 1)   # mL

amount_water_aragonite = calculate_volume_pure("water", 0.21)   # mL

wash_water = calculate_volume_pure("water", 0.5)   # mL
n_washes_water = 6
wash_water_once = wash_water/n_washes_water

wash_acetone = 0.3*wash_water   # mL
n_washes_acetone = 3
wash_acetone_once = wash_acetone/n_washes_acetone

reaction_volume = amount_CaCl2 + amount_K2CO3

# ===============================================================================

#universal constants
default_speed = 40   # mL/min
default_speed_wash = 70   # mL/min
dead_volume_filter = 10   # mL, see 'prime_filter'

#reaction-specific constants
reactor_stirring_rate = 300   # RPM
reaction_temp = 30   # 'C



###########################
#    Defined  Functions   #
###########################

def wash_left (flask_washing_solvent:str, wash_volume:float=3, repeats:int=1):
    '''Washes backbone left of flask_washing_solvent'''
    for _ in range(repeats): 
        c.move(flask_washing_solvent, "waste1", wash_volume, speed=default_speed_wash, use_backbone=True)

def wash_right (flask_washing_solvent:str, wash_volume:float=3, repeats:int=1):
    '''Washes backbone right of flask_washing_solvent'''
    for _ in range(repeats): 
        c.move(flask_washing_solvent, "waste3", wash_volume, speed=default_speed_wash, use_backbone=True)

def wash_backbone (flask_washing_solvent:str, wash_volume:float=3):
    '''Washes backbone in its entirety with the solvent in an indicated flask'''
    wash_left(flask_washing_solvent, wash_volume)
    wash_right(flask_washing_solvent, wash_volume)
    print("wash_backbone complete")

def prime_inputs(inputs_to_be_primed, with_volume):
    '''Used to prime inputs before synthesis to reduce error in initial transfer'''
    for i in inputs_to_be_primed:
        c.move(i, "waste1", with_volume)

def prime_filter(solvent_bottom:str) -> None:
    '''Primes filter for use by filling up to the frit with solvent from indicated input'''
    c.move(solvent_bottom, "filter", dead_volume_filter, speed=default_speed, dest_port="bottom")
    c.wait(15)   # secs, pause for 15 secs
    print("prime_filter complete")

def vac_dry (time):
    # turn to vacuum
    c['valve1'].move_to_position(5)
    c.wait(time)
    c['valve1'].move_to_position(0)



###########################
#        Procedure        #
###########################

# PRIOR TO PROCEDURE START:
linebreak = "\n***************************************************************************\n"
print(linebreak, "PREPARATION OF CALCIUM CARBONATE POLYMORPHS", linebreak)

# choose polymorph and raise error if unacceptable value passed
polymorph = input("Please select polymorph to prepare: (1: Vaterite, 2: Aragonite, 3: Calcite)")

if polymorph in ['1','2','3']:
    pass
else:
    raise ValueError("Invalid Entry. Valid values are '1', '2', '3' not '%s'." % polymorph)

# ask for required reagents
print("\nSOLUTIONS/SOLVENTS REQUIRED:")
print("aqueous calcium chloride dihydrate (input12, %.2f g, %.2fM, %.2f mL)" % (calculate_mass("CaCl2.2H2O", 1), float(reagent_conc["CaCl2.2H2O"]), amount_CaCl2))
print("aqueous potassium carbonate (input14, %.2f g, %.2fM, %.2f mL)" % (calculate_mass("K2CO3", 1), float(reagent_conc["K2CO3"]), amount_K2CO3))
if polymorph=='1':
    print("water (input20, %.2f mL)" % wash_water)
    print("acetone (input19, %.2f mL)" % wash_acetone)
elif polymorph=='2':
    print("water (input20, %.2f mL - with %.2f mL for reflux)" % ((amount_water_aragonite + wash_water + wash_water), amount_water_aragonite))
    print("acetone (input19, %.2f mL)" % (wash_acetone + wash_acetone))
elif polymorph=='3':
    print("water (input20, %.2f mL)" % wash_water)
    print("acetone (input19, %.2f mL)" % wash_acetone)
c.breakpoint()

# setup camera
c.start_recording(0)
c.camera.change_recording_speed(20)

# prime inputs
prime_inputs(["input20", "input12", "input20", "input19", "input20", "input14", "input20", "input19"], 3)
wash_backbone("input20", 3)



# Ca1 VATERITE SYNTHESIS
###########################

# prepare filter
prime_filter("input20")
c.move("input14", "filter", amount_K2CO3, dest_port="top")
c.move("input5", "filter", 5, dest_port="top")
wash_backbone("input20")
wash_left("input20")
wash_left("input20")
c.stirrer.set_stir_rate("filter", 300)   # rpm
c['julabo_chiller'].set_temperature(reaction_temp)   # °C
c.stirrer.stir("filter")
c['julabo_chiller'].start() 
c.wait(900)   # secs, 15 mins, time to allow contents of filter to reach correct temperature


# addition of CaCl2.2H2O
c.move("input12", "filter", amount_CaCl2, dest_port="top", end_pump_speed=25)
c.move("input5", "filter", 5, dest_port="top", end_pump_speed=25)


# prepare for waiting period and sediment
wash_backbone("input20")
wash_backbone("input20")
c.stirrer.stop_stir("filter")
c['julabo_chiller'].stop()


# sedimentation, filtration, and washing w/water and acetone
if polymorph in ['1','2']:
    c['julabo_chiller'].set_temperature(reaction_temp)   # °C
    c['julabo_chiller'].start() 
    c.wait(420)   # secs, 7 mins
    
    c.move("filter", "waste1", (2*reaction_volume), src_port="bottom")

    for i in range(n_washes_water):
        c.move("input20", "filter", wash_water_once)
        c.wait(30)
        c.move("filter", "waste1", wash_water_once*2, initial_pump_speed=25)
        print("Completed water wash #", i+1, " of ", n_washes_water, ".")

    vac_dry(5)

    for i in range(n_washes_acetone):
        c.move("input19", "filter", wash_acetone_once)
        c.wait(30)
        c.move("filter", "waste1", wash_acetone_once*2, initial_pump_speed=25)
        print("Completed acetone wash #", i+1, " of ", n_washes_acetone, ".")

    vac_dry(5)

    c['julabo_chiller'].stop()

    c.move("input5", "filter", 5)
    c.move("filter", "waste1", 30, initial_pump_speed=25)

if polymorph in ['1']:
    print("\nCa1 Vaterite: reaction end.", linebreak)



# Ca2 ARAGONITE SYNTHESIS
###########################
if polymorph in ['2']:
    # reaction setup and reflux
    wash_left("input20")
    wash_left("input20")
    c.move("input20", "filter", amount_water_aragonite)

    c.stirrer.set_stir_rate("filter", 300)   # rpm
    c['julabo_chiller'].set_temperature(100)   # °C
    c.stirrer.stir("filter")
    c['julabo_chiller'].start() 
    c.wait(7200)   # secs, 2 hrs
    c['julabo_chiller'].stop()
    c.stirrer.stop_stir("filter")


    # sedimentation, filtration, and washing w/water and acetone
    c.wait(600)   # secs, 10 mins
    
    c.move("filter", "waste1", (2*amount_water_aragonite), src_port="bottom")

    for i in range(n_washes_water):
        c.move("input20", "filter", wash_water_once)
        c.wait(30)
        c.move("filter", "waste1", wash_water_once*2, initial_pump_speed=25)
        print("Completed water wash #", i+1, " of ", n_washes_water, ".")

    for i in range(n_washes_acetone):
        c.move("input19", "filter", wash_acetone_once)
        c.wait(30)
        c.move("filter", "waste1", wash_acetone_once*2, initial_pump_speed=25)
        print("Completed acetone wash #", i+1, " of ", n_washes_acetone, ".")

    c.move("input5", "filter", 5)
    c.move("filter", "waste1", 30, initial_pump_speed=25)

    print("\nCa2 Aragonite: reaction end.", linebreak)



# Ca3 CALCITE SYNTHESIS
###########################
if polymorph in ['3']:
    # reflux in mother liquor
    c.stirrer.set_stir_rate("filter", 300)   # rpm
    c['julabo_chiller'].set_temperature(60)   # °C
    c.stirrer.stir("filter")
    c['julabo_chiller'].start() 
    c.wait(18000)   # secs, 5 hrs
    c['julabo_chiller'].stop()
    c.stirrer.stop_stir("filter")


    # sedimentation, filtration, and washing w/water and acetone
    #c.breakpoint()   # need to define sedimentation time - I expect the breakpoint to fail
    c.wait(7200)
    
    c.move("filter", "waste1", (2*reaction_volume), src_port="bottom")

    for i in range(n_washes_water):
        c.move("input20", "filter", wash_water_once)
        c.wait(30)
        c.move("filter", "waste1", wash_water_once*2, initial_pump_speed=25)
        print("Completed water wash #", i+1, " of ", n_washes_water, ".")

    for i in range(n_washes_acetone):
        c.move("input19", "filter", wash_acetone_once)
        c.wait(30)
        c.move("filter", "waste1", wash_acetone_once*2, initial_pump_speed=25)
        print("Completed acetone wash #", i+1, " of ", n_washes_acetone, ".")

    c.move("input5", "filter", 5)
    c.move("filter", "waste1", 30, initial_pump_speed=25)

    print("\nCa3 Calcite: reaction end.", linebreak)



# CLEANING
###########################
# clean backbone
wash_backbone("input20")
wash_backbone("input20")
wash_backbone("input19")



###########################
#        Reset All        #
###########################

c.camera.stop()

# resets working directory
os.chdir("C:/Users/group/")
print("Reset working directory to: " + os.getcwd())

