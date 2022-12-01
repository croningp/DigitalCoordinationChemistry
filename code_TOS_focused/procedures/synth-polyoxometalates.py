##########################################################################################################
# ACC_RuLibrary1 - Reaction
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
GRAPH_FILE = path + "\\graph-general.json"

# create Chempiler object
c = Chempiler(
    experiment_code = "polyoxometalate", 
    output_dir=path, 
    simulation=True, 
    graph_file=GRAPH_FILE, 
    device_modules=[ChemputerAPI, commanduinolabware]
)



###########################
#   User-Set Variables    #
###########################

limiting_reagent_moles = 0.001608   # mol, limiting reagent: Na2MoO4.2H2O
average_wash = 10   # mL, used for washes after crystallisation



###########################
#   Automatic Variables   #
###########################

#define basic calculations for amounts of substance
def calculate_mass(reagent, equiv):
    MW = reagent_MW.get(str(reagent))
    mass = float(MW) * (float(limiting_reagent_moles) * float(equiv))
    return round(mass, 2)

def calculate_volume_soln(reagent, equiv):
    conc = reagent_conc.get(str(reagent))
    volume_L = (float(limiting_reagent_moles) * float(equiv)) / float(conc)
    volume_mL = volume_L * 1000
    return round(volume_mL, 2)

def calculate_volume_pure(reagent, conc_wrt_LR):   
    volume_L = limiting_reagent_moles / conc_wrt_LR   # L
    volume_mL = volume_L * 1000
    return round(volume_mL, 2)

#define dictionaries containing physical data on reagents
reagent_MW = {
    "Na2MoO4.2H2O":"241.95",
    "NaOAc":"82.03",
    "AcOH":"60.05",
    "water":"18.02", 
    "Na2S2O4":"174.11",
}   # g/mol
reagent_conc = {
    "NaOAc":"0.975",
    "AcOH":"8.7",   # 50% when glacial acetic acid is 17.4 M
    "HCl":"1"
}   # Molar
reagent_density = {
    "water":"1.000",
}   # g/mL

#calculate amounts of substances
amount_molybdate = calculate_mass("Na2MoO4.2H2O", 1)   # g
amount_dithionite = calculate_mass("Na2S2O4", 0.35)
amount_NaOAc = calculate_volume_soln("NaOAc", 6.06)   # mL
amount_AcOH = calculate_volume_soln("AcOH", 27)   # mL
amount_Brown_HCl = 0   # mL
amount_Blue_HCl = calculate_volume_soln("HCl", 0)   # mL
amount_water = calculate_volume_pure("water", 0.3216)   # mL

reaction_volume = amount_NaOAc + amount_AcOH + amount_Blue_HCl + amount_water

# ===============================================================================

#universal constants
default_speed = 40   # mL/min
default_speed_wash = 70   # mL/min
dead_volume_filter = 10   # mL, see 'prime_filter'

#reaction-specific constants
reactor_stirring_rate = 300   # RPM



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

def clean_reactors(reactors:list, solvent="input20", wash_volume=10, temp=100) -> None:
    '''
    Cleans each reactor ready for next set of reactions. 
    '''
    c.camera.change_recording_speed(100)

    for i in reactors:
        c.move(solvent, i, wash_volume)

    for i in reactors:
        c.stirrer.set_stir_rate(i, 300)
        c.stirrer.stir(i)
        c.stirrer.set_temp(i, temp)
        c.stirrer.heat(i)
        
    c.wait(600)   # secs, 10 mins
    
    for i in reactors:
        c.stirrer.stop_heat(i)
    
    for i in reactors:
        c.move(i, "waste3", wash_volume*2)
    
    for i in reactors:
        c.move(solvent, i, wash_volume)
    
    c.wait(60)   # secs, 1 min

    for i in reactors:
        c.move(i, "waste3", wash_volume*2)
    
    for i in reactors:
        c.move(solvent, i, wash_volume)
    
    c.wait(60)   # secs, 1 min

    for i in reactors:
        c.move(i, "waste3", wash_volume*2)

    for i in reactors:
        c.stirrer.stop_stir(i)



###########################
#        Procedure        #
###########################

# PRIOR TO PROCEDURE START:
linebreak = "\n***************************************************************************\n"
print(linebreak, "PREPARATION OF MOLYBDENUM BLUES AND BROWNS", linebreak)

print("\nSOLID REAGENTS REQUIRED:")
print("sodium molybdate dihydrate (r1, %.2f g)" % amount_molybdate)
print("sodium dithionite (r2, %.2f g)" % amount_dithionite)
print("\nSOLUTIONS/SOLVENTS REQUIRED:")
print("aqueous sodium acetate (input9, %.2f mL, %.2fM)" % (amount_NaOAc, float(reagent_conc["NaOAc"])))
print("50%% acetic acid (input10, %.2f mL, %.2fM)" % (amount_AcOH, float(reagent_conc["AcOH"])))
print("For Mo Blue only: 1M HCl (input11, %.2f mL, %.2fM)" % (amount_Blue_HCl, float(reagent_conc["HCl"])))
print("water (input20, %.2f mL)" % amount_water)
print("acetone (input19)")
c.breakpoint()

c.start_recording(0)

# prime inputs
c.camera.change_recording_speed(20)
prime_inputs(["input9", "input10", "input11"], 3)
wash_backbone("input20", 3)


# prepare reaction mix in r1
c.move("input9", "r1", 10, )   # mL, 0.975M sodium acetate(aq)
c.move("input10", "r1", 5)   # mL, 50% acetic acid(aq)
c.move("input11", "r1", 0)    # mL, 1M HCl(aq)
c.move("input20", "r1", 5)   # mL, water
c.move("input5", "r1", 5)   # air to blow in XS

wash_backbone("input20")

c.stirrer.set_stir_rate("r1", 300)   # rpm
c.stirrer.stir("r1")
c.wait(30)   # secs
c.stirrer.stop_stir("r1")


# combine with solid sodium dithionite in r2
c.move("r1", "r2", 40)   # mL
c.move("input5", "r2", 5)   # air to blow in XS

wash_backbone("input20")


# reaction
c.stirrer.set_stir_rate("r2", 300)   # rpm
c.stirrer.stir("r2")
c.wait(600)   # secs
c.stirrer.stop_stir("r2")


# sample reaction
c.move("r2", "uvsample", 50, end_pump_speed=15)
c.move("input5", "uvsample", 5)   # clear tubing with air


# clean rig
wash_backbone("input20")
wash_backbone("input20")
wash_backbone("input20")

c.move("input5", "uvsample", 10, speed=60)
clean_reactors(["r1","r2","r3","r4"], solvent="input20", wash_volume=10, temp=50)
wash_backbone("input20")


###########################
#        Reset All        #
###########################

c.camera.stop()

# resets working directory
os.chdir("C:/Users/group/")
print("Reset working directory to: " + os.getcwd())
