import matplotlib.pyplot as plt
import csv


############################################################## ACCESSING AME2020 DATA FROM CSV ########################################################################

mass_numbers = []                                                          # Holds each nucleus’s mass number
binding_energies = []                                                      # Holds binding energy per nucleon
estimated_flags = []                                                       # Records whether the CSV file marked that row as estimated

with open("nuclear_data.csv", newline="", encoding="utf-8") as file:       # Opening the CSV file

    reader = csv.DictReader(file)                                          # Reads each row using column names
    
    for row in reader:

        mass_numbers.append(int(row["A"]))                                 # Saves that row's mass number
        binding_energies.append(float(row["BE_per_nucleon_MeV"]))          # Saves binding energy per nucleon as a decimal number
        estimated_flags.append(row["estimated"] == "True")                 # Saves True/False as seen in CSV file

not_estimated_A = []                                                       # Stores mass numbers for rows not flagged as estimated 
not_estimated_BE = []                                                      # Stores binding energies per nucleon for rows not flagged as estimated
estimated_A = []                                                           # Stores mass numbers for rows flagged as estimated
estimated_BE = []                                                          # Stores binding energies per nucleon for rows flagged as estimated


for A, BE, estimated in zip(mass_numbers, binding_energies, estimated_flags):            # Pairs each corresponding A, B/A, and flag

    if estimated:
        estimated_A.append(A)                                                            # Adding to the list, A for CSV estimated flag true
        estimated_BE.append(BE)                                                          # Adding to the list, BE/A for CSV estimated flag true

    else:
        not_estimated_A.append(A)                                                        # Adding to the list, A for CSV estimated flag false
        not_estimated_BE.append(BE)                                                      # Adding to the list, BE/A for CSV estimated flag false


highest_binding_energy_per_nucleon_by_A= {}                                              # Creates an empty dictionary

for A, BE in zip(mass_numbers, binding_energies):

    if A not in highest_binding_energy_per_nucleon_by_A or BE > highest_binding_energy_per_nucleon_by_A[A]:

        highest_binding_energy_per_nucleon_by_A[A]= BE                                   # Finds or replaces values for highest binding energy per nucleon for each

curve_A= sorted(highest_binding_energy_per_nucleon_by_A)
curve_BE= [highest_binding_energy_per_nucleon_by_A[A] for A in curve_A]


########################################################################## DEFINING SEMF ###############################################################################

a_v= 15.75                                        # Volume coefficient
a_s= 17.8                                         # Surface coefficient
a_c= 0.711                                        # Coulomb coefficient
a_a= 23.7                                         # Asymmetry coefficient
a_p= 34.0                                        # Pairing coefficient




semf_A= []                                        # Empty list for mass number
semf_BE= []                                       # Empty list for binding energy per nucleon


for A in range(16, max(mass_numbers) + 1):  

    best_model_value = 0.0                        # Starts the search for the largest predicted BE/A for each A

    for Z in range(1, A):

        N = A - Z                                 # Calculating neutron number

        pairing = 0.0                             # Starting at zero; also even-odd condition

        if Z % 2 == 0 and N % 2 == 0:             # Even condition

            pairing = a_p / A**(3 / 4)

        elif Z % 2 == 1 and N % 2 == 1:           # Odd condition

            pairing = -a_p / A**(3 / 4)

        if A == 56 and Z == 26:
            print("Even-even check:", pairing)
        if A == 56 and Z == 25:
            print("Odd-odd check:", pairing)
        if A == 55 and Z == 26:
            print("Odd-A check:", pairing)

        total_BE = ( a_v * A - a_s * A**(2 / 3) - a_c * Z * (Z - 1) / A**(1 / 3) - a_a * (A - 2 * Z)**2 / A + pairing )       # SEMF equation

        model_value = total_BE / A                # Converts total binding energy into binding energy per nucleon

        if model_value > best_model_value:
            best_model_value = model_value        # Keeps the largest predicted value as we try different Z values

    semf_A.append(A)                              # Adds the mass number
    semf_BE.append(best_model_value)              # Adds or replaces binding energy per nucleon




difference_A = []                                                                   # Stores how much difference is in SEMF and AME2020
difference_BE = []

for A, model_BE in zip(semf_A, semf_BE):                                           
    if A in highest_binding_energy_per_nucleon_by_A:
        data_BE = highest_binding_energy_per_nucleon_by_A[A]                        # Calculates the difference in MeV per nucleon. 
        difference_A.append(A)                                                      # A positive result means the data’s upper envelope is above the SEMF prediction
        difference_BE.append(data_BE - model_BE)                                    # A negative result means it is below the SEMF prediction




########################################################################## PLOTTING ####################################################################################

plt.figure(figsize= (12,6))
plt.title("Nuclear binding energy curve")
plt.xlabel("Mass number, A")
plt.ylabel("Binding energy per nucleon (MeV)")

plt.scatter(not_estimated_A, not_estimated_BE, s=8, alpha=0.22, label="Not marked estimated")             # Plots nuclei for which CSV estimated flag is false
plt.scatter(estimated_A, estimated_BE, s=8, alpha=0.22, label="Marked estimated")                         # Plots nuclei for which CSV estimated flag is true

plt.plot(curve_A, curve_BE, color="black", linewidth=1.5, label="Upper envelope (all CSV entries)")       # Adds a black line connecting BE/A for each A

plt.plot(semf_A, semf_BE, color="red", linestyle="--", linewidth=3, zorder=5, label="SEMF with pairing, A >= 16")      # Adds red dashed line for SEMF prediction

plt.legend()
plt.grid(True)                                                                           # Draws grids


plt.savefig("Binding_energy_curve(CSV data).png")                                        # Binding energy curve from AME2020 data
plt.show()




