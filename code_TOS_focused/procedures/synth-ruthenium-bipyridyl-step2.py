##########################################################################################################
# Ruthenium bipyridyl (Step 2)
##########################################################################################################
from chempiler import Chempiler
import ChemputerAPI
import appdirs
import os

# set working directory to same as this script
os.chdir(os.path.dirname(os.path.abspath(__file__)))
path=os.getcwd()

# define graph
GRAPH_FILE = path + "\\general-graph.json"
XDL_FILE = 0


# Create Chempiler object
c = Chempiler(
    experiment_code = "ruthenium-bipyridyl-2", 
    output_dir=path, 
    simulation=True, 
    graph_file=GRAPH_FILE, 
    device_modules=[ChemputerAPI]
)

###########################
#   User-Set Variables    #
###########################

limiting_reagent_moles = 0.0005   # mol, limiting reagent: Ru(bipy)2Cl2.xH2O
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
    "RuCl3.xH2O":"261.47", 
    "LiCl":"42.39",
    "bipyridyl":"156.188", 
    "DMF":"73.095", 
    "acetone":"58.080",
    "Ru(bipy)2Cl2.2H2O":"520.38",   # assumed to be dihydrate (also exists as de- or monohydrate)
    "EtOH":"46.069"
    }   # g/mol
reagent_conc = {
    "LiCl":"3",   # in DMF
    "bipyridyl1":"1.5",   # in DMF
    "bipyridyl2":"0.2"   # in EtOH
    }   # Molar
reagent_density = {
    "DMF":"0.948", 
    "acetone":"0.7845",
    "EtOH":"0.7893"
    }   # g/mL

#calculate amounts of substances
amount_Rubipy2Cl2 = calculate_mass("Ru(bipy)2Cl2.2H2O", equiv=1)
amount_bipyridyl = calculate_volume_soln("bipyridyl2", equiv=1)
amount_EtOH = round(calculate_volume_pure("EtOH", conc_wrt_LR=0.03) - amount_bipyridyl)

amount_EtOH_wash = average_wash   # mL, see 'User-Set Variables'
amount_Et2O_wash = average_wash   # mL, see 'User-Set Variables'

reaction_volume = amount_bipyridyl + amount_EtOH
amount_aq_sat_NaBF4 = 0.08 * reaction_volume

# ===============================================================================

#universal constants
default_speed = 40   # mL/min
default_speed_wash = 70   # mL/min
dead_volume_filter = 10   # mL, see 'prime_filter'

#reaction-specific constants
reactor_stirring_rate = 150   # RPM
dissolution_temp = 75   # °C
reaction_temp = 85   # °C
cryst_temp = 4   # °C


###########################
#    Defined  Functions   #
###########################

def prime_filter(solvent_bottom):
    # primes filter for use by filling up to the frit with solvent from indicated flask
    c.move(solvent_bottom, "filter", dead_volume_filter, speed=default_speed, dest_port="bottom")
    c.wait(15)   # secs, pause for 15 secs
    print("prime_filter complete")

def wash_precipitate(flask_washing_solvent, wash_volume, filter_to_wash="filter", destination="waste1", stir_speed=60):
    # washes residue from filtration that remains on the frit with a specified volume of solvent.
    # also requires specification of the waste the solvent is pushed to
    c.move(flask_washing_solvent, filter_to_wash, wash_volume, speed=default_speed_wash, dest_port="top")
    if wash_volume%6 !=0:
        n_pulls = wash_volume//6 + 1
    else:
        n_pulls = wash_volume//6
    if isinstance(n_pulls, int) == True:
        c.move(filter_to_wash, destination, volume=((n_pulls+1)*10), speed=20, src_port="bottom")
        c.wait(15)
        c.move(filter_to_wash, destination, volume=10, speed=20, src_port="bottom")
        print("wash_precipitate complete")

def prime_inputs(inputs_to_be_primed, with_volume):
    # used to prime inputs before synthesis
    for i in inputs_to_be_primed:
        c.move(i, "waste1", with_volume)

def wash_left (flask_washing_solvent, wash_volume=3, repeats=1):
    # Washes backbone left of flask_washing_solvent
    for _ in range(repeats): 
        c.move(flask_washing_solvent, "waste1", wash_volume, speed=default_speed_wash, use_backbone=True)

def wash_right (flask_washing_solvent, wash_volume=3, repeats=1):
    # Washes backbone right of flask_washing_solvent
    for _ in range(repeats): 
        c.move(flask_washing_solvent, "waste3", wash_volume, speed=default_speed_wash, use_backbone=True)

def wash_backbone (flask_washing_solvent, wash_volume=3):
    # Washes backbone in its entirety with the solvent in an indicated flask
    wash_left(flask_washing_solvent, wash_volume)
    wash_right(flask_washing_solvent, wash_volume)
    print("wash_backbone complete")


###########################
#        Procedure        #
###########################

# PRIOR TO PROCEDURE START:
linebreak = "\n***************************************************************************\n"
print(linebreak, "PREPARATION OF TRIS(BIPYRIDINE)RUTHENIUM(II) TETRAFLUOROBORATE", linebreak)
print("SOLID REAGENTS REQUIRED:")
print("bis(2,2'-bipyridine)dichlororuthenium(II) dihydrate (filter, %.2f g)" % calculate_mass("Ru(bipy)2Cl2.2H2O", equiv=1))
print("\nSOLUTIONS/SOLVENTS REQUIRED:")
print("2,2-bipyridyl soln in EtOH (input14, %s M, %.2f mL, %.2f g)"  % (reagent_conc.get("bipyridyl2"), amount_bipyridyl, calculate_mass("bipyridyl", equiv=1.00))) 
print("EtOH (input13, 5x %.2f mL)" % amount_EtOH)
print("aqueous saturated sodium tetrafluoroborate (input12, %.2f mL)" % amount_aq_sat_NaBF4)
print("water (input20)")
print("diethyl ether (input18)")
print("\nCLEANING SOLVENTS REQUIRED: aqueous 0.1M HCl (input21), aqueous 0.05M Na4EDTA (input22), acetonitrile (input23)")
c.breakpoint()

# start recording
c.start_recording(0)
c.camera.change_recording_speed(20)

# prime backbone and inputs
prime_inputs(["input14","input12","input20","input19","input18"], 3)   # don't need to do with input13, as backbone wash achieves this
wash_backbone("input13", 3)   # w/ EtOH

# dissolve complex
prime_filter("input13")
c.move("input13", "filter", amount_EtOH-2)   # EtOH (note that 2 mL is held back to wash bipyridyl into flask later on)
c['julabo_chiller'].set_temperature(reaction_temp)   # °C
c.stirrer.set_stir_rate("filter", reactor_stirring_rate)
c.stirrer.stir("filter")
c['julabo_chiller'].start()
c.camera.change_recording_speed(100)
c.wait(900)   # secs, 15 mins   gives some time for dissolution before the addition of other reagents
#c.stirrer.wait_for_temp("filter")

# reaction
c.camera.change_recording_speed(20)
c.move("input14", "filter", amount_bipyridyl)   # 2,2'-bipyridyl soln in EtOH
c.move("input13", "filter", 2)   # wash in to filter with held over EtOH
wash_backbone("input13", 3)   # w/EtOH
c.camera.change_recording_speed(100)
for i in range(4):
    c.wait(17280)   # secs, 4.8 hrs - repeated 4 times for a total of 19.2 hrs
    c.move("input13", "filter", amount_EtOH)
c.wait(17280)   # secs, 4.8 hrs - added to the above gives total rxn time of 24 hrs
c['julabo_chiller'].stop()
c.stirrer.stop_stir("filter")

# precipitation and filtration
c.camera.change_recording_speed(20)
prime_inputs(["input12"], 3)
c.stirrer.set_stir_rate("filter", 30)
c.stirrer.stir("filter")   # stir gently with overhead stirrer to mix properly
c.move("input12", "filter", amount_aq_sat_NaBF4, dest_port="top")   # sat NaBF4(aq)
c.move("input5", "filter", 3, dest_port="top")  # use air to blow in the remaining sat NaBF4(aq)
wash_backbone("input13", 3)   # w/EtOH
c.wait(600)   # secs, 10 mins to permit precipitation
c.stirrer.stop_stir("filter")

# sanity check
print("Begin Recryst?")
c.breakpoint()

# recrystallisation
prime_filter("input13")
c.move("input13", "filter", (100 - 10 - amount_aq_sat_NaBF4 - 2 * reaction_volume))
c.stirrer.set_stir_rate("filter", 150)
c['julabo_chiller'].set_temperature(100)
c.stirrer.stir("filter")
c['julabo_chiller'].start()
c.breakpoint()   # to permit heating for as long as is required to dissolve all material
c.stirrer.stop_stir("filter")
c['julabo_chiller'].stop()

c.breakpoint()   # to permit precipitation to occur

# filtration and wash residue
c.move("filter", "waste1", 100, src_port="bottom")   # filter
wash_precipitate("input13", amount_EtOH_wash)   # w/ EtOH
wash_precipitate("input13", amount_EtOH_wash)   # w/ EtOH
wash_precipitate("input18", amount_Et2O_wash)   # w/ diethyl ether

# clean backbone
wash_left("input20", repeats=5)
c.camera.stop()


###########################
#        Reset All        #
###########################

# resets stirrer-hotplate settings
for i in ["r1","r2","r3","r4"]:
    c.stirrer.set_stir_rate(i, 0)
    c.stirrer.stop_stir(i)
    c.stirrer.set_temp(i, 0)
    c.stirrer.stop_heat(i)

# resets overhead stirrer
c.stirrer.set_stir_rate("filter", 0)
c.stirrer.stop_stir("filter")

# resets chiller
c['julabo_chiller'].set_temperature(20)
c['julabo_chiller'].stop()

# resets working directory
os.chdir("C:/Users/group/")
print("Reset working directory to: " + os.getcwd())