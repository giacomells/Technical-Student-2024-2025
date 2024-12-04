#DICTONARY OF THE VARIABLES FOR THE NOTEBOOK

p = 400.0  # beam momentum (GeV/c)
momentum = p  # beam momentum (GeV/c)
Brho = p * 3.3356  # beam rigidity ???

N_EX = 10e-6
N_EY = 5e-6
DPP = 1e-4
# Set the septum aperture size
septum_aperture_position = 68e-3  # Size of the septum aperture in meters

        
# Names of the elements
septum_namesZS = ["zs.21633", "zs.21639", "zs.21655", "zs.21671", "zs.21676"]
septum_namesMST = ["mst.21774", "mst.21779", "mst.21794"]
septum_namesMSE = ["mse.21832", "mse.21837", "mse.21852", "mse.21857", "mse.21872"]



# ZS Features
ZSentrances = [0., 3.91, 7.82, 11.73, 15.64]
ZSexits = [3.1, 7.01, 10.92,  14.83, 18.74]
blade_positionsZS = [68e-3, 62.12e-3, 56.23e-3, 50.35e-3, 44.46e-3]
heightZSelements = 4.665e-3 # Height in x
widthZSelements = 3.1  # Width of the entire element in x    
apertureZS = 2.0e-2  #Aperture between the blade and the cathode 
kickZS = 0.440e-3 / 5

# TPST Features
blade_positionTPST = 40e-3
thicknessTPST = 5.2e-3
heightTPST = 1.06e-3
widthTPST = 2.14
thicknessTPST = 4.2e-3
apertureTPST = 40e-3
kickTPST = 0

# MST Features
MSTentrances = [2.79, 6.02, 9.25]
MSTexits = [5.19, 8.42, 11.65]
blade_positionsMST = [40.79e-3, 42.40e-3, 44.01e-3]
heightMSTelements = 1.19e-3 # height in x
widthMSTelements =  2.4 # Width of the entire element in x
thicknessMST = 4.2e-3
apertureMST = 97.9e-3
kickMST = 1.69520713e-3 / 3

#  MSE Features
blade_positionsMSEentrance = [52.70e-3, 51.84e-3, 55.11e-3, 62.50e-3, 74.00e-3]
blade_positionsMSEexits = [51.67e-3, 53.72e-3, 59.89e-3, 70.19e-3, 84.60e-3]
thicknessMSE = 17.25e-3
widthMSEelements = 2.28
apertureMSE = 61e-3
kickMSE = 9.74519477e-3 / 5