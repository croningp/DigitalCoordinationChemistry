##########################################################################################################
# Library Generation via Late-Stage Variation
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
GRAPH_FILE = path + "\\graph-library-generation.json" 

simulation_state=False

# create Chempiler object
c = Chempiler(
    experiment_code = "library-generation", 
    output_dir=path, 
    simulation=simulation_state, 
    graph_file=GRAPH_FILE, 
    device_modules=[ChemputerAPI, commanduinolabware]
)

# logging
#logging.basicConfig(level=logging.DEBUG)

# import sample wheel
cfg = {
    "ios": [
        {
            "type": "serial",
            "port": "COM6"
        },
    ],
    "devices": {
        "wheel": {
            "command_id": "X",
            "config": {
                "reverted_switch": True,
                "reverted_direction": False,
                "enabled_acceleration": False,
                "speed": 5000,
                "homing_speed": 12000,
                "acceleration": 2000
            }
        }
    }
}

if simulation_state == True:
    # code for simulated modular wheel device
    device=c.graph.graph.nodes['modular_wheel']['obj']
else:
    # code for real modular wheel device
    c.graph.graph.nodes['modular_wheel']['obj'] = commanduinolabware.CommanduinoLabware(config=cfg)



###########################
#  Import .csv Rxn List   #
###########################

# import input list of reactions to be carried out and randomise
rxn_list = pd.read_csv("library-generation.csv", index_col=0)
rxn_list_rand = rxn_list

# create list of reactors the same length as the number of reactions to be performed
reactor_list = ["r1","r2","r3","r4"]
df_Reactor = []
n, i = 0, 0

while n <= (len(rxn_list_rand)-1):
    next_addition = reactor_list[i]
    df_Reactor.append(next_addition)
    if i == (len(reactor_list)-1):
        i = 0
    else:
        i += 1
    n += 1

# combine reactor assignment with randomised reaction list
rxn_list_rand['Reactor']=df_Reactor

# save randomised reaction list as .csv
now = datetime.datetime.now()
file_name = "library-generation_date%2.0f%2.0f_time%.0f%.0f.csv" % (now.month, now.day, now.hour, now.minute)
rxn_list_rand.to_csv(path_or_buf=file_name)



###########################
#  Set-up Email Updates   #
###########################

port = 465  # For SSL
smtp_server = "smtp.gmail.com"
sender_email = ""  # Enter your address
receiver_email = ""  # Enter receiver address
message = """\
Subject: Library Generation Experiment Update

A set of reactions has completed: """

def email_update(df, password):
    '''
    Sends email update to 'receiver_email' defined above 
    indicating that the set of reactions recorded in 'df' have completed.
    Uses a specially set-up gmail account as an smtp server.
    '''
    ls = []

    for index, row in df.iterrows():
        variable = ("Rxn " + str(index) + " in " + row['Reactor'])
        ls.append(variable)
        
    context = ssl.create_default_context()

    with smtplib.SMTP_SSL(smtp_server, port, context=context) as server:
        server.login(sender_email, password)
        server.sendmail(sender_email, receiver_email, (message+str(ls)))



###########################
#   Automatic Variables   #
###########################

# universal constants
default_speed = 40   # mL/min
default_speed_wash = 70   # mL/min
dead_volume_filter = 10   # mL, see 'prime_filter'

# ===============================================================================

linebreak = "\n***************************************************************************\n"



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

def test_hardware(reactors:list, filterstir:bool=True, filterheat:bool=True) -> None:
    '''
    Triggers reactors, modular wheel, and filter stirrer and/or heater individually.
    Acts as confirmation that each is functioning correctly.
    '''
    c.move("input5", "waste3", 3, speed=70)

    for i in reactors:
        c.stirrer.set_stir_rate(i, 300)   # rpm
        c.stirrer.stir(i)
        c.stirrer.set_temp(i, 50)   # °C
        c.stirrer.heat(i)
    c.wait(30)
    for i in reactors:
        c.stirrer.stop_stir(i)
        c.stirrer.stop_heat(i)

    c['modular_wheel'].turn_motor("wheel", 1)
    c['modular_wheel'].turn_motor("wheel", 10)

    if filterstir == True:
        c.stirrer.set_stir_rate("filter", 150)   # rpm
        c.stirrer.stir("filter")
        c.wait(30)
        c.stirrer.stop_stir("filter")

    if filterheat == True:
        c['julabo_chiller'].set_temperature(75)   # °C
        c['julabo_chiller'].start()
        c.wait(30)
        c['julabo_chiller'].stop()

    print("Test Completed. Continue with experiment?")
    c.breakpoint()

def dissolve_Ru1(volume:float=60, duration:int=600) -> None:
    '''
    Dissolves Ru(bpy)2Cl2 intermediate in specified volume of ethanol
    Heats to 75 °C for specified duration to ensure full dissolution
    '''
    c.move("input13", "filter", volume)   # mL, EtOH
    c['julabo_chiller'].set_temperature(75)   # °C
    c.stirrer.set_stir_rate("filter", 150)   # rpm
    c.stirrer.stir("filter")
    c['julabo_chiller'].start()
    c.camera.change_recording_speed(100)
    c.wait(duration)    # secs, 10 mins   gives some time for dissolution before the addition of other reagents
    c.stirrer.stop_stir("filter")
    c['julabo_chiller'].stop()

def reaction_setup(df, reactors:list) -> None:
    '''
    Transfers Ru1 to each reactor in use from filter.
    Transfers ligands to each reactor as requested (via .csv file).
    '''
    print(linebreak, '\n', df)

    c.camera.change_recording_speed(20)

    # ensure metal complex solution is homogeneous
    c.stirrer.set_stir_rate("filter", 150)
    c.stirrer.stir("filter")
    c.wait(20)   # secs, 1 min

    for i in reactors:
        c.move("filter", i, 4)   # transfer Ru1
        c.move("input5", i, 5)   # blow in with air

    c.stirrer.stop_stir("filter")
    wash_backbone("input13", 2)
    wash_backbone("input13", 2)
    wash_backbone("input13", 2)

    # adds components specified in the .csv file
    for index, row in df.iterrows():
        c.move(str(row['Component1']), str(row['Reactor']), 4)   # add first component
        c.move("input5", row['Reactor'], 5)   # blow in with air
        c.move(str(row['Component2']), str(row['Reactor']), 4)   # add second component
        c.move("input5", str(row['Reactor']), 5)   # blow in with air
        wash_backbone("input13", 2)
        wash_backbone("input13", 2)

def run_reactions(reactors:list) -> None:
    '''
    Set and start each reactor heating and stirring.
    Switch off again after for a 16 hour period.
    '''
    c.camera.change_recording_speed(100)

    for i in reactors:
        c.stirrer.set_stir_rate(i, 300)   # rpm
        c.stirrer.stir(i)
        c.stirrer.set_temp(i, 85)   # °C
        c.stirrer.heat(i)

    c.wait(57600)   # secs, 16 hrs

    # stop heating and stirring
    for i in reactors:
        c.stirrer.stop_stir(i)
        c.stirrer.stop_heat(i)

def sample_reactions(reactors:list) -> None:
    '''
    Move product mix to the modular wheel.
    Turn wheel ready to receive next sample.
    '''
    c.camera.change_recording_speed(20)

    for i in reactors:
        c.move(i, "modular_wheel", 20)
        c.move("input13", "modular_wheel", 3)   # wash in with EtOH
        c.move("input5", "modular_wheel", 5)   # blow in with air
        c['modular_wheel'].turn_motor("wheel", 1)
        wash_backbone("input13")

def clean_reactors(reactors:list, wash_volume=10) -> None:
    '''
    Cleans each reactor ready for next set of reactions. 
    '''
    c.camera.change_recording_speed(100)

    for i in reactors:
        c.move("input13", i, wash_volume)

    for i in reactors:
        c.stirrer.set_stir_rate(i, 300)
        c.stirrer.stir(i)
        c.stirrer.set_temp(i, 85)
        c.stirrer.heat(i)
        
    c.wait(600)   # secs, 10 mins
    
    for i in reactors:
        c.stirrer.stop_heat(i)
    
    for i in reactors:
        c.move(i, "waste3", wash_volume*2)
    
    for i in reactors:
        c.move("input13", i, wash_volume)
    
    c.wait(60)   # secs, 1 min

    for i in reactors:
        c.move(i, "waste3", wash_volume*2)
    
    for i in reactors:
        c.move("input13", i, wash_volume)
    
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
print(linebreak, "GENERATION OF A LIBRARY OF BIS(2,2'-BIPYRIDINE) RUTHENIUM PYRIDYLS", linebreak)

#password = getpass.getpass("Please provide password for automated smtp server:")

print("\nSOLID REAGENTS REQUIRED:")
print("Bis(2,2'-bipyridine) ruthenium(II) chloride (filter, ~1 g !!Ensure weight recorded!!)")
print("\nSOLUTIONS/SOLVENTS REQUIRED:")
print("6x 0.06M solns of pyridyl derivatives in ethanol: 2,6-lutidine (input8), "
"2-aminopyridine (input9), 2-picolinic acid (input10), isoquinoline (input11), "
"piperazine (input12), and 4,4'-bipyridine (input14)")
print("ethanol (input13)")

print("\nRANDOMISED LIST OF REACTIONS: \n", rxn_list_rand)


# hardware test
print("\n", "On continuing, hardware will activate to test, and then will prompt operator to continue.")
c.breakpoint()
test_hardware(reactor_list, filterstir=False, filterheat=False)

c.start_recording(0)

# prime inputs
c.camera.change_recording_speed(20)
prime_inputs(["input8", "input9", "input10", "input11", "input12", "input14"],3)
wash_backbone("input13", 3)

# dissolve Ru1 intermediate
prime_filter("input13")
dissolve_Ru1(volume=60, duration=300)

# remove excess solvent used to prime filter
c.camera.change_recording_speed(20)
c.move("filter", "waste1", 10)   

# reactions 1-4
rxns2 = rxn_list_rand.iloc[0:4]
### during run: dissolve_Ru1(0, 300)
reaction_setup(rxns2, reactor_list)
run_reactions(reactor_list)
sample_reactions(reactor_list)
clean_reactors(reactor_list)
#email_update(rxns2, password)

# reactions 5-8
rxns3 = rxn_list_rand.iloc[4:8]
### during run: dissolve_Ru1(0, 300)
reaction_setup(rxns3, reactor_list)
run_reactions(reactor_list)
sample_reactions(reactor_list)
clean_reactors(reactor_list)
#email_update(rxns3, password)


###########################
#        Reset All        #
###########################
c.camera.stop()

# resets working directory
os.chdir("C:/Users/group/")
print("Reset working directory to: " + os.getcwd())
