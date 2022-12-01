##########################################################################################################
# Cisplatin (Step 2)
##########################################################################################################
from chempiler import Chempiler
import ChemputerAPI
import appdirs
import os

# set working directory to same as this script
os.chdir(os.path.dirname(os.path.abspath(__file__)))
path=os.getcwd()

# define graph
GRAPH_FILE = path + "general-graph.json"
XDL_FILE = 0


# Create Chempiler object
c = Chempiler(
    experiment_code = "cisplatin-2", 
    output_dir=path, 
    simulation=True, 
    graph_file=GRAPH_FILE, 
    device_modules=[ChemputerAPI]
)


###########################
#   User-Set Variables    #
###########################

limiting_reagent_moles = 0.000870   # mol, limiting reagent: PtI2(NH3)2
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
    "PtI2(NH3)2":"482.95",
    "AgNO3":"169.87",
    "water":"18.02", 
    "KCl":"74.55",
}   # g/mol
reagent_conc = {}   # Molar
reagent_density = {
    "water":"1.000",
}   # g/mL

#calculate amounts of substances
amount_cisdiiododiammineplatinum = calculate_mass("PtI2(NH3)2", 1)   # g
amount_AgNO3 = calculate_mass("AgNO3", 2)   # g
amount_KCl = calculate_mass("KCl", 3)   # g
amount_water = calculate_volume_pure("water", 0.25)   # mL

amount_water_wash = average_wash   # mL, see 'User-Set Variables'
amount_EtOH_wash = average_wash   # mL, see 'User-Set Variables'
amount_Et2O_wash = average_wash   # mL, see 'User-Set Variables'

reaction_volume = amount_water

# ===============================================================================

#universal constants
default_speed = 40   # mL/min
default_speed_wash = 70   # mL/min
dead_volume_filter = 10   # mL, see 'prime_filter'

#reaction-specific constants
reactor_stirring_rate = 150   # RPM
filter_stirring_rate = 150   # RPM
reaction_temp = 40   # °C
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
print(linebreak, "PREPARATION OF CISPLATIN, [Pt(NH3)2Cl2]", linebreak)
print("SOLID REAGENTS REQUIRED:")
print("cis-diamminediiodoplatinum(II) (r1, %.2f g)" % amount_cisdiiododiammineplatinum)
print("silver nitrate (r1, %.2f g)" % amount_AgNO3)
print("potassium chloride (filter, %.2f g)" % amount_KCl)
print("\nSOLUTIONS/SOLVENTS REQUIRED:")
print("water (input20, %.2f mL)"  % amount_water) 
print("ethanol (input13)")
print("diethyl ether (input18)")
print("additional water (input20)")
print("\nCLEANING SOLVENTS REQUIRED: acetone (input19), aqueous 0.1M HCl (input21), aqueous 0.05M Na4EDTA (input22), acetonitrile (input23)")
c.breakpoint()

# start recording
c.start_recording(0)
c.camera.change_recording_speed(20)

# prime backbone and inputs
wash_backbone("input20", 3)
prime_inputs(["input19", "input18", "input13"], 3)
wash_backbone("input20", 3)

# reaction 1 - substitution of iodide ligands for nitrates
c.move("input20", "r1", amount_water)
c.move("input5", "r1", 5, speed=5)   # air to blow in water, slowly so it doesn't spray the reaction mix everywhere
c.stirrer.set_stir_rate("r1", reactor_stirring_rate)
c.stirrer.stir("r1")
c.wait(57600)   # secs, 16 hrs
c.stirrer.stop_stir("r1")
c.breakpoint()

# prime filter
c.move("input20", "filter", dead_volume_filter, dest_port="bottom")

# filter using 'in-line' filter (with reactor washing)
c.stirrer.stir("r1")
c.move("input20", "r1", 5)
c.move("r1", "filter", amount_water + 10, speed=5)
c.move("input20", "r1", 10)
c.move("r1", "filter", 15, speed=5)
c.move("input20", "r1", 10)
c.move("r1", "filter", 15, speed=5)
check = input("Flush Again? y/N")
check2 = input("With water? y/N")
while check == "y":
    if check2 == "y":
        c.move("input20", "r1", 10)
    c.move("r1", "filter", 15, speed=5)
    check = input("Wash Again? y/N")
    check2 = input("With water? y/N")
c.stirrer.stop_stir("r1")
wash_right("input20")
wash_right("input20")

# reaction 2 - substitution of nitrates for chloride ligands
c.stirrer.set_stir_rate("filter", filter_stirring_rate)
c.stirrer.stir("filter")
c['julabo_chiller'].set_temperature(reaction_temp)   # °C
c['julabo_chiller'].start() 
c.wait(900)   # secs, 15 mins
c['julabo_chiller'].stop()
c.stirrer.stop_stir("filter")
c['julabo_chiller'].set_temperature(cryst_temp)
c['julabo_chiller'].start()
c.wait(900)   # secs, 15 mins
c['julabo_chiller'].stop()

c.breakpoint()   # wait until permitted to continue

# filter and wash
c.move("filter", "cryst", volume=90, speed=20, src_port="bottom")
wash_precipitate("input13", amount_EtOH_wash)
wash_precipitate("input18", amount_Et2O_wash)

# clean backbone
wash_backbone("input20")
wash_backbone("input20")
wash_backbone("input19")
c.camera.stop()



###########################
#        Reset All        #
###########################

# resets working directory
os.chdir("C:/Users/group/")
print("Reset working directory to: " + os.getcwd())
