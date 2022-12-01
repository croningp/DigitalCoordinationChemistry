##########################################################################################################
# Basic Manganese Acetate, {Mn3O}
##########################################################################################################
from chempiler import Chempiler
import ChemputerAPI
import appdirs
import os

# set working directory to same as this script
os.chdir(os.path.dirname(os.path.abspath(__file__)))
path=os.getcwd()

# define graph
GRAPH_FILE = path+"general-graph.json"
XDL_FILE = 0

# Create Chempiler object
c = Chempiler(
    experiment_code = "basic-manganese-acetate", 
    output_dir=path, 
    simulation=True, 
    graph_file=GRAPH_FILE, 
    device_modules=[ChemputerAPI]
)

###########################
#   User-Set Variables    #
###########################

limiting_reagent_moles = 0.00407   # mol, limiting reagent: Mn(OAc)2.4H2O
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
    "Mn(OAc)2.4H2O":"245.09",
    "Bu4NMnO4":"361.4",
    "EtOH":"46.07", 
    "pyridine":"79.10",
    "AcOH":"60.05", 
    "NaClO4":"122.44" 
}   # g/mol
reagent_conc = {
    "NaClO4":"0.95"   # approx for saturated NaClO4 given solubility of 14.7g/100g EtOH @ 25'C
                      # Source: http://chemister.ru/Database/properties-en.php?dbid=1&id=769
}   # Molar
reagent_density = {
    "EtOH":"0.789",
    "pyridine":"0.978",
    "AcOH":"1.049"
}   # g/mL

#calculate amounts of substances
amount_Bu4NMnO4 = calculate_mass("Bu4NMnO4", equiv=0.405)
amount_solvent = calculate_volume_pure("", conc_wrt_LR=0.23)   # volume of entire solvent mixture
amount_EtOH = round((amount_solvent/35)*20)   # 20:3:12 EtOH/pyr/AcOH
amount_pyr = round((amount_solvent/35)*3)
amount_AcOH = round((amount_solvent/35)*12)
amount_satNaClO4 = calculate_volume_soln("NaClO4", equiv=0.7)

amount_EtOH_wash = average_wash   # mL, see 'User-Set Variables'
amount_Et2O_wash = average_wash   # mL, see 'User-Set Variables'
amount_acetone = 30   # mL, see 'User-Set Variables'

reaction_volume = amount_EtOH + amount_pyr + amount_AcOH

# ===============================================================================

#universal constants
default_speed = 40   # mL/min
default_speed_wash = 70   # mL/min
dead_volume_filter = 10   # mL, see 'prime_filter'

#reaction-specific constants
reactor_stirring_rate = 300   # RPM
filter_stirring_rate = 150   # RPM
cryst_temp = 60   # °C


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
print(linebreak, "PREPARATION OF [Mn3O(OAc)6(pyr)3](ClO4)", linebreak)
print("SOLID REAGENTS REQUIRED:")
print("manganese(II) acetate tetrahydrate (r1, %.2f g)" % calculate_mass("Mn(OAc)2.4H2O", equiv=1))
print("n-butylammonium permanganate (filter, %.2f g)" % amount_Bu4NMnO4)
print("\nSOLUTIONS/SOLVENTS REQUIRED:")
print("EtOH (input13, %.2f mL)"  % amount_EtOH) 
print("pyridine (input9, %.2f mL)" % amount_pyr)
print("glacial acetic acid (input10, %.2f mL)" % amount_AcOH)
print("saturated ethanolic NaClO4 (input6, %s M, %.2f mL, %.2f g)" % (reagent_conc.get("NaClO4"), amount_satNaClO4, calculate_mass("NaClO4", equiv=0.7)))
print("water (input20)")
print("acetone (input19)")
print("diethyl ether (input18)")
print("\nCLEANING SOLVENTS REQUIRED: aqueous 0.1M HCl (input21), aqueous 0.05M Na4EDTA (input22), acetonitrile (input23)")
c.breakpoint()

# start recording
c.start_recording(0)
c.camera.change_recording_speed(20)

# prime backbone and inputs
prime_inputs(["input9","input10","input18","input20", "input19"], 3)   # don't need to do with input13, as backbone wash achieves this
wash_backbone("input13", 3)   # w/ EtOH

# dissolve complex
c.move("input13", "r1", (amount_EtOH))   # EtOH (2 mL held back to wash)
c.move("input9", "r1", amount_pyr)   # pyridine
c.move("input10", "r1", amount_AcOH)   # AcOH
c.move("input5", "r1", 5)   # blow in with air
c.stirrer.set_stir_rate("r1", reactor_stirring_rate)
c.stirrer.stir("r1")
c.wait(600)   # secs, 10 mins
c.stirrer.stop_stir("r1")

# combination with butylammonium permanganate
c.move("r1", "filter", (reaction_volume+10))
prime_filter("input13")   # primed after combining w/ rxn mix to avoid decomposition
c.stirrer.set_stir_rate("filter", filter_stirring_rate)
c.stirrer.stir("filter")
c.move("r1", "filter", 5)   # blow in any remaining rxn mix
wash_backbone("input13", 3)   # w/EtOH
c.wait(900)   # secs, 15 mins

# add sat NaClO4 soln
prime_inputs(["input6"], 3)
c.move("input6", "filter", amount_satNaClO4)   # sat NacClO4(aq)
c.move("input5", "filter", 5)   # blow in with air
wash_backbone("input13", 3)   # w/EtOH
c.wait(600)   # secs, 10 mins
c.stirrer.stop_stir("filter")

# precipitation
c.wait(10800)   # secs, 3 hrs

# filtration and wash residue
c.move("filter", "waste1", reaction_volume + amount_satNaClO4 + dead_volume_filter + 3, src_port="bottom")   # filter
wash_precipitate("input13", amount_EtOH_wash)   # w/ EtOH
wash_precipitate("input18", amount_Et2O_wash)   # w/ diethyl ether
wash_precipitate("input18", amount_Et2O_wash)   # w/ diethyl ether

# recrystallisation
prime_filter("input19")
c.move("input19", "filter", amount_acetone)
c.stirrer.set_stir_rate("filter", filter_stirring_rate)
c['julabo_chiller'].set_temperature(cryst_temp)   # °C
c.stirrer.stir("filter")
c['julabo_chiller'].start()
c.wait(1200)   # secs, 20 mins
c.stirrer.stop_stir("filter")
c.wait(7200)   # secs, 2 hrs
c['julabo_chiller'].stop()
c['julabo_chiller'].set_temperature(4)
c['julabo_chiller'].start()
c.wait(3600)   # secs, 1 hrs
c['julabo_chiller'].stop()

# filtration and washing
c.move("filter", "cryst", (amount_acetone+30), src_port="bottom", speed=20)
wash_precipitate("input18", amount_Et2O_wash)   # w/ diethyl ether

# clean backbone
wash_backbone("input20")
wash_backbone("input20")
wash_backbone("input20")
c.camera.stop()


###########################
#        Reset All        #
###########################

# resets working directory
os.chdir("C:/Users/group/")
print("Reset working directory to: " + os.getcwd())