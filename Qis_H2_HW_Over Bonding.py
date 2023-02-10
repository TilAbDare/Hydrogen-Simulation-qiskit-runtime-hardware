# https://qiskit.org/documentation/nature/tutorials/03_ground_state_solvers.html
# --------------------------------------------------------------------
# ******************  Importing libraries ****************************
# --------------------------------------------------------------------
import timeit
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from qiskit import IBMQ
from qiskit.algorithms.optimizers import SPSA
from qiskit.circuit.library import EfficientSU2
from qiskit_nature.runtime import VQEClient
from qiskit_nature.second_q.algorithms import GroundStateEigensolver
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import ParityMapper
from qiskit_nature.second_q.mappers import QubitConverter
from qiskit_nature.second_q.transformers import ActiveSpaceTransformer
from qiskit_nature.units import DistanceUnit

start = timeit.default_timer()
np.set_printoptions(precision=4, suppress=True)

# --------------------------------------------------------------------
# **************************  IBMQ Access ****************************
# --------------------------------------------------------------------

#IBMQ.save_account('ba9ad96cbaecbb73dc94d234b5ced5a2422de52d088705cac4676f2e606c99b4179af4583758f611f6392fe21d9e179d8180a2013b605016556013b2dbf50797')
#IBMQ.load_account()
#IBMQ.providers()
#provider = IBMQ.get_provider(hub="ibm-q")  # replace by your runtime provider
#backend = provider.get_backend("ibmq_qasm_simulator")  # select a backend that supports the runtime

api_key = None
for attempt in range(3):                                                # <- (1)
    try:                                                                # <- (3)
        if api_key:
            IBMQ.save_account("ba9ad96cbaecbb73dc94d234b5ced5a2422de52d088705cac4676f2e606c99b4179af4583758f611f6392fe21d9e179d8180a2013b605016556013b2dbf50797", overwrite=True)
        else:
            IBMQ.load_account()
        break                                                           # <- (4)
    except Exception:
        api_key = input("Enter IBMQ API Key (attempt %s): " % attempt)  # <- (2)


#IBMQ.providers()
#provider = IBMQ.get_provider(hub="ibm-q")  # replace by your runtime provider
#backend = provider.get_backend("ibmq_qasm_simulator")  # select a backend that supports the runtime


# --------------------------------------------------------------------
# ************************* Molecule Setup ***************************
# --------------------------------------------------------------------
h2_length = np.arange(0.2, 2.9, 0.05)
hw_energy_list = []

for bond_length in h2_length:
    driver = PySCFDriver(
        atom='H 0.0 0.0 0.0; H 0.0 0.0 '+str(bond_length),
        charge=0,
        spin=0,
        unit=DistanceUnit.ANGSTROM,
        basis='sto3g',
    )
    problem = driver.run()
    active_space_trafo = ActiveSpaceTransformer(
        num_electrons=problem.num_particles, num_spatial_orbitals=2)
    # transform the electronic structure problem
    problem = active_space_trafo.transform(problem)
    qubit_converter = QubitConverter(ParityMapper(), two_qubit_reduction=True)
    ansatz = EfficientSU2(num_qubits=2, reps=1, entanglement="linear", insert_barriers=True)
    optimizer = SPSA(maxiter=100)
    initial_point = np.random.random(ansatz.num_parameters)
    runtime_vqe = VQEClient(
        ansatz=ansatz,
        optimizer=optimizer,
        initial_point=initial_point,
        provider=provider,
        backend=backend,
        shots=1,
        measurement_error_mitigation=True,
    )
    runtime_vqe_groundstate_solver = GroundStateEigensolver(qubit_converter, runtime_vqe)
    runtime_vqe_result = runtime_vqe_groundstate_solver.solve(problem)
    hw_energy_list += [runtime_vqe_result.total_energies[0]]

#ansatz.decompose().draw("mpl", style="iqx")


stop = timeit.default_timer()
runtime = stop - start
print('Run Time: ', runtime, 'sec', ' or ', runtime/60, 'min','\n')
# --------------------------------------------------------------------
# ********************* Saving The Data locally **********************
# --------------------------------------------------------------------
df = pd.DataFrame(list(zip(h2_length,  hw_energy_list)),
                  columns = ['distance', 'IBM Hardware'])
#df.to_excel("Diag.xlsx")
print('-'*8, ' Final Results ', '-'*8)
print(df)

# --------------------------------------------------------------------
# ******************************  Plot *******************************
# --------------------------------------------------------------------

plt.plot(h2_length, hw_energy_list, 'o', color='red', markerfacecolor=(1, 1, 0, 0.5), label='IBM Hardware', markersize = 7)
plt.title("Potential Energy Curve of hydrogen Molecule")
plt.xlabel("H-H band length")
plt.ylabel("Energy Hartree")
plt.legend()
#plt.show()


# --------------------------------------------------------------------
# ********************* Saving The Data locally **********************
# --------------------------------------------------------------------

df = pd.DataFrame(list(zip(h2_length, hw_energy_list)),
                  columns = ['h2_length', 'hw_energy_list'])
df.to_excel("Qis_H2_Hardware.xlsx")
print(df)




















