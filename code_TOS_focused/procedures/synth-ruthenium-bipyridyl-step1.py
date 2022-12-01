##########################################################################################################
# Ruthenium bipyridyl (Step 1)
##########################################################################################################

from chempiler import Chempiler
import ChemputerAPI
import appdirs
import os

# set working directory to same as this script
os.chdir(os.path.dirname(os.path.abspath(__file__)))
path=os.getcwd()

# define graph
GRAPH_FILE = path + "\\graph-general.json"
XDL_FILE = 0

# create Chempiler object
c = Chempiler(
    experiment_code = "ruthenium-pyridyl-1", 
    output_dir=appdirs.user_data_dir("xdl"), 
    simulation=True, 
    graph_file=GRAPH_FILE, 
    device_modules=[ChemputerAPI]
)

###########################
#   User-Set Variables    #
###########################

limiting_reagent_moles = 0.003   # mol, limiting reagent: RuCl3.xH2O (lit. 0.00182)
average_wash = 20   # mL, used for washes after crystallisation

###########################
#   Automatic Variables   #
###########################

# define basic calculations for amounts of substance
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

# define dictionaries containing physical data on reagents
reagent_MW = {
    "RuCl3.xH2O":"261.47", 
    "LiCl":"42.39", 
    "bipyridyl":"156.188", 
    "DMF":"73.095", 
    "acetone":"58.080"
    }   # g/mol
reagent_conc = {
    "LiCl":"2", 
    "bipyridyl":"1.1"
    }   # Molar
reagent_density = {
    "DMF":"0.948", 
    "acetone":"0.7845"
    }   # g/mL

# calculate amounts of substances
amount_RuCl3 = calculate_mass("RuCl3.xH2O", equiv=1)
amount_bipyridyl = calculate_volume_soln("bipyridyl", equiv=2.00)
amount_LiCl = calculate_volume_soln("LiCl", equiv=7.25)
amount_DMF = round(calculate_volume_pure("DMF", conc_wrt_LR=0.182) - amount_bipyridyl - amount_LiCl, 2)
amount_acetone = calculate_volume_pure("acetone", conc_wrt_LR=0.061)
amount_water_wash = average_wash   # mL, see 'User-Set Variables'
amount_Et2O_wash = average_wash   # mL, see 'User-Set Variables'

reaction_volume = amount_bipyridyl + amount_LiCl + amount_DMF

# ===============================================================================

# universal constants
default_speed = 40   # mL/min
default_speed_wash = 70   # mL/min
dead_volume_filter = 10   # mL, see 'prime_filter'

# reaction-specific constants
reactor_stirring_rate = 150   # RPM
reaction_temp = 160   # °C
cryst_temp = 4   # °C


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

def wash_backbone (flask_washing_solvent, wash_volume):
    # Washes backbone in its entirety with the solvent in an indicated flask
    c.move(flask_washing_solvent, "waste1", wash_volume, speed=default_speed_wash)
    c.move(flask_washing_solvent, "waste3", wash_volume, speed=default_speed_wash)
    print("wash_backbone complete")

def wash_reactor (flask_washing_solvent, reactor_to_wash, wash_volume, destination):
    # washes an individual reactor with the solvent in an indicated flask
    # also requires specification of the waste the solvent is pushed to
    c.move(flask_washing_solvent, reactor_to_wash, wash_volume, speed=default_speed_wash)
    c.stirrer.set_stir_rate(reactor_to_wash, 600)
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

###########################
#        Procedure        #
###########################

# PRIOR TO PROCEDURE START:
print("SOLID REAGENTS REQUIRED:")
print("ruthenium trichloride hydrate (filter, %.2f g)" % calculate_mass("RuCl3.xH2O", equiv=1))
print()
print("SOLUTIONS/SOLVENTS REQUIRED:")
print("2,2-bipyridyl soln in DMF (input15, %s M, %.2f mL, %.2f g)"  % (reagent_conc.get("bipyridyl"), amount_bipyridyl, calculate_mass("bipyridyl", equiv=2.00))) 
print("lithium chloride soln in DMF (input16, %s M, %.2f mL, %.2f g)"  % (reagent_conc.get("LiCl"), amount_LiCl, calculate_mass("LiCl", equiv=7.25)))
print("DMF (input17, %.2f mL)" % amount_DMF)
print("acetone (input19, %.2f mL)" % amount_acetone)
print("water (input20)")
print("diethyl ether (input18)")
c.breakpoint()

# start recording
c.start_recording(0)
c.camera.change_recording_speed(20)

# prime backbone and inputs
prime_inputs(["input15","input16","input20","input19", "input18"], 3)   # don't need to do with input17, as backbone wash achieves this
wash_backbone("input17", 3)   # w/ DMF

# set-up
prime_filter("input17")   # w/ DMF
c.move("input15", "filter", amount_bipyridyl, dest_port="top")   # 2,2-bipyridyl soln in DMF
c.move("input16", "filter", amount_LiCl, dest_port="top")   # LiCl soln in DMF
c.move("input17", "filter", amount_DMF, dest_port="top")   # DMF
c.move("input5", "filter", 5, dest_port="top")   # blow in with air
wash_backbone("input17", 3)   # w/ DMF

# reaction
c.stirrer.set_stir_rate("filter", reactor_stirring_rate)
c['julabo_chiller'].set_temperature(reaction_temp)   # °C
c.stirrer.stir("filter")
c['julabo_chiller'].start()
c.camera.change_recording_speed(100)
c.wait(10800)   # secs, 3 hrs
c['julabo_chiller'].stop()
c.stirrer.stop_stir("filter")

# cool to r.t.
c['julabo_chiller'].set_temperature(30)
c.wait(10800)    # secs, 3 hrs
#c.stirrer.wait_for_temp("filter")    # can I use wait_for_temp with c['julabo_chiller']?

# crystallisation
c.camera.change_recording_speed(20)
c['julabo_chiller'].set_temperature(cryst_temp)   # °C
c.stirrer.set_stir_rate("filter", 30)
c['julabo_chiller'].start()
c.stirrer.stir("filter")
c.move("input19", "filter", amount_acetone, dest_port="top")   # acetone
wash_backbone("input19", 3)   # w/ acetone
c.camera.change_recording_speed(100)
c.wait(32400)   # secs, 9 hrs
c.camera.change_recording_speed(25)
c['julabo_chiller'].stop()
c.stirrer.stop_stir("filter")

# filtration
c.move("filter", "waste1", reaction_volume + dead_volume_filter + amount_acetone + 3, src_port="bottom")

# washing residue
wash_precipitate("input20", amount_water_wash, "waste1")   # w/ water
wash_precipitate("input20", amount_water_wash, "waste1")   # w/ water
wash_precipitate("input18", amount_water_wash, "waste1")   # w/ diethyl ether

wash_backbone("input19", 3)
wash_backbone("input20", 3)
wash_backbone("input18", 3)

# finish - prompt to collect crystals and wash rig
print("Process finished. Please collect crystals. Wash rig? (not including filter)")
c.breakpoint()

print("Process finished. Wash filter?")
c.breakpoint()

c.camera.stop()


###########################
#        Reset All        #
###########################

# resets working directory
os.chdir("C:/Users/group/")
print("Reset working directory to: " + os.getcwd())