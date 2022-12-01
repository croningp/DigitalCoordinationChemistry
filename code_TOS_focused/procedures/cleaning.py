##########################################################################################################
# ACC_General_Cleaning
##########################################################################################################
from chempiler import Chempiler
import ChemputerAPI
import appdirs
import os

# set working directory to same as this script
os.chdir(os.path.dirname(os.path.abspath(__file__)))
path=os.getcwd()

# define graph
GRAPH_FILE = path + "\graph-cleaning.json"
XDL_FILE = 0


# Create Chempiler object
c = Chempiler(
    experiment_code = "Maintenance_Cleaning", 
    output_dir=appdirs.user_data_dir("xdl"), 
    simulation=True, 
    graph_file=GRAPH_FILE, 
    device_modules=[ChemputerAPI]
)

###########################
#   Automatic Variables   #
###########################

#universal constants
default_speed = 40   # mL/min
default_speed_wash = 60   # mL/min
dead_volume_filter = 10   # mL, see 'prime_filter'

###########################
#    Defined  Functions   #
###########################

def prime_filter(solvent_bottom):
    # primes filter for use by filling up to the frit with solvent from indicated flask
    c.move(solvent_bottom, "filter", dead_volume_filter, speed=default_speed, dest_port="bottom")
    c.wait(15)   # secs, pause for 15 secs

def wash_precipitate(flask_washing_solvent, volume_washing_solvent, destination):
    # washes residue from filtration that remains on the frit with a specified volume of solvent.
    # also requires specification of the waste the solvent is pushed to
    c.move(flask_washing_solvent, "filter", volume_washing_solvent, speed=default_speed)
    c.wait(60)   # secs, soaking for 1 min
    c.move ("filter", destination, (1.5 * volume_washing_solvent), speed=default_speed, src_port="bottom")
    c.wait(15)   # secs, pause for 15 secs

def prime_inputs(inputs_to_be_primed, with_volume):
    # used to prime inputs before synthesis
    for i in inputs_to_be_primed:
        c.move(i, "waste1", with_volume)

def wash_left (flask_washing_solvent, wash_volume=3, repeats=1):
    # Washes backbone left of flask_washing_solvent
    for i in range(repeats): 
        c.move(flask_washing_solvent, "waste1", wash_volume, speed=default_speed_wash, use_backbone=True)

def wash_right (flask_washing_solvent, wash_volume=3, repeats=1):
    # Washes backbone right of flask_washing_solvent
    for i in range(repeats): 
        c.move(flask_washing_solvent, "waste3", wash_volume, speed=default_speed_wash, use_backbone=True)

def wash_backbone (flask_washing_solvent, wash_volume=3):
    # Washes backbone in its entirety with the solvent in an indicated flask
    wash_left(flask_washing_solvent, wash_volume)
    wash_right(flask_washing_solvent, wash_volume)
    print("wash_backbone complete")

def wash_reactor (flask_washing_solvent, reactor_to_wash, destination, wash_volume=10, stirrer_speed=400):
    # washes an individual reactor with the solvent in an indicated flask
    # also requires specification of the waste the solvent is pushed to
    c.move(flask_washing_solvent, reactor_to_wash, wash_volume, speed=default_speed_wash)
    c.stirrer.set_stir_rate(reactor_to_wash, stirrer_speed)
    c.stirrer.stir(reactor_to_wash)
    c.wait(60)   # secs
    c.stirrer.stop_stir(reactor_to_wash)
    c.move(reactor_to_wash, destination, (wash_volume + 10), speed=default_speed_wash)
    print("wash_reactor complete: ", reactor_to_wash)
    
def wash_filter (flask_washing_solvent, filter_to_wash, wash_volume, destination):
    # washes the filter (first from top and then from bottom port) with solvent in an indicated flask
    # also requires specification of the waste the solvent is pushed to
    c.move(flask_washing_solvent, filter_to_wash, wash_volume, speed=default_speed_wash, dest_port="top")
    c.wait(45)   # secs
    c.move(filter_to_wash, destination, (wash_volume + 10), speed=default_speed_wash, src_port="bottom")
    print("wash_filter complete: ", filter_to_wash)

def wash_reactor_hot (flask_washing_solvent, reactor_to_wash, destination, wash_volume=10, stirrer_speed=400, at_temp=20):
    # washes an individual reactor with the solvent in an indicated flask by refluxing at given temperature
    # also requires specification of the waste the solvent is pushed to
    c.move(flask_washing_solvent, reactor_to_wash, wash_volume, speed=default_speed_wash)
    c.stirrer.set_stir_rate(reactor_to_wash, stirrer_speed)
    c.stirrer.set_temp(reactor_to_wash, at_temp)
    c.stirrer.heat(reactor_to_wash)
    c.stirrer.stir(reactor_to_wash)
    c.wait(600)   # secs
    c.stirrer.stop_stir(reactor_to_wash)
    c.stirrer.stop_heat(reactor_to_wash)
    c.move(reactor_to_wash, destination, (wash_volume + 10), speed=default_speed_wash)
    print("wash_reactor_hot complete: ", reactor_to_wash)

###########################
#        Procedure        #
###########################

print("Wash?")
c.breakpoint()

wash_reactor("input20", "r1", "waste2")   # water

wash_reactor_hot("input19", "r1", "waste2", at_temp=65)   # reflux w/acetone
wash_reactor("input19", "r1", "waste2")   # acetone
wash_reactor("input21", "r1", "waste2")   # 0.1M HCl(aq): 2.1 mL fuming made up to 250 mL
wash_reactor("input22", "r1", "waste2")   # 0.05M EDTA4-(aq): 4.75 g made up to 250 mL
wash_right("input20", 3, 2)   # water
wash_backbone("input19")   # acetone

wash_left("input21", repeats=3)
wash_left("input22", repeats=1)
wash_backbone("input20")
wash_backbone("input19")
