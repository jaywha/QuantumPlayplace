# Variational Quantum Algorithms - Cirq Tutorial
# Followed on July 27 2018
# C&P by Jay Whaley
#
# Function for error printing to avoid clogging the output from SO
# Spinning loading "bar" to keep track of long-operations from GitHub (Halo)

from __future__ import print_function
from halo import Halo
import cirq
import random
import sys
import numpy as np


# Standard Error Stream Printing: from https://stackoverflow.com/a/14981125/7476183
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)
    pass


pass

# Prepare Quantum Circuit for Operations
# Also shows moments, time-sliced rotations, and grid placement strategies.

grid_length = 3  # define the length of the grid.
qubits = [cirq.GridQubit(i, j) for i in range(grid_length) for j in range(grid_length)]  # define qubits on the grid.
eprint(qubits)
eprint("\n")

circuit = cirq.Circuit()
circuit.append(cirq.H.on(q) for q in qubits if (q.row + q.col) % 2 == 0)
circuit.append(cirq.X(q) for q in qubits if (q.row + q.col) % 2 == 1)
eprint(circuit)
eprint("\n")

for i, m in enumerate(circuit):
    eprint('Moment {}: {}'.format(i, m))
    pass
eprint("\n")

circuit = cirq.Circuit()
circuit.append([cirq.H.on(q) for q in qubits if (q.row + q.col) % 2 == 0],
               strategy=cirq.InsertStrategy.EARLIEST)
circuit.append([cirq.X(q) for q in qubits if (q.row + q.col) % 2 == 1],
               strategy=cirq.InsertStrategy.EARLIEST)
eprint(circuit)
eprint("\n")


# Define Relevant Rotations

def rot_x_layer(length, half_turns):
    """Yields X rotations by half_turns on a square grid of given length."""
    rot = cirq.RotXGate(half_turns=half_turns)
    for i in range(length):
        for j in range(length):
            yield rot(cirq.GridQubit(i, j))
            pass
        pass
    pass


pass


def rot_z_layer(field_h, half_turns):
    """Yields Z rotations by half_turns conditioned on the field h."""
    gate = cirq.RotZGate(half_turns=half_turns)
    for i, h_row in enumerate(field_h):
        for j, h_ij in enumerate(h_row):
            if h_ij == 1:
                yield gate(cirq.GridQubit(i, j))
                pass
            pass
        pass
    pass


pass


def rot_11_layer(el_jr, el_jc, half_turns):
    """Yields rotations about |11> conditioned on the jr and jc fields."""
    gate = cirq.Rot11Gate(half_turns=half_turns)
    for i, jr_row in enumerate(el_jr):
        for j, jr_ij in enumerate(jr_row):
            if jr_ij == -1:
                yield cirq.X(cirq.GridQubit(i, j))
                yield cirq.X(cirq.GridQubit(i + 1, j))
            yield gate(cirq.GridQubit(i, j),
                       cirq.GridQubit(i + 1, j))
            if jr_ij == -1:
                yield cirq.X(cirq.GridQubit(i, j))
                yield cirq.X(cirq.GridQubit(i + 1, j))
                pass
            pass
        pass
    pass

    for i, jc_row in enumerate(el_jc):
        for j, jc_ij in enumerate(jc_row):
            if jc_ij == 1:
                yield gate(cirq.GridQubit(i, j),
                           cirq.GridQubit(i, j + 1))
                pass
            pass
        pass
    pass


pass  # rot_11_layer

circuit = cirq.Circuit()
circuit.append(rot_x_layer(2, 0.1))
eprint(circuit)
eprint("\n")


# Generate random problem instances

def rand2d(rows, cols):
    return [[random.choice([+1, -1]) for _ in range(rows)] for _ in range(cols)]


pass


def random_instance(length):
    # transverse field terms
    ri_h = rand2d(length, length)
    # links within a row
    ri_jr = rand2d(length, length - 1)
    # links within a column
    ri_jc = rand2d(length - 1, length)
    return ri_h, ri_jr, ri_jc


pass


def one_step(o_h, o_jr, o_jc, x_half_turns, h_half_turns, j_half_turns):
    length = len(o_h)
    yield rot_x_layer(length, x_half_turns)
    yield rot_z_layer(o_h, h_half_turns)
    yield rot_11_layer(o_jr, o_jc, j_half_turns)
    pass


pass


# Calculate Energy
def energy_func(length, e_h, e_jr, e_jc):
    def energy(measurements):
        meas_list_of_lists = [measurements[i:i + length] for i in range(length)]
        pm_meas = 1 - 2 * np.array(meas_list_of_lists).astype(np.int32)
        tot_energy = np.sum(pm_meas * e_h)
        for i, jr_row in enumerate(e_jr):
            for j, jr_ij in enumerate(jr_row):
                tot_energy += jr_ij * pm_meas[i, j] * pm_meas[i + 1, j]
        for i, jc_row in enumerate(e_jc):
            for j, jc_ij in enumerate(jc_row):
                tot_energy += jc_ij * pm_meas[i, j] * pm_meas[i, j + 1]
        return tot_energy

    return energy


pass


# Objective Function
def obj_func(obj_result):
    energy_hist = obj_result.histogram(key='x', fold_func=energy_func(3, h, jr, jc))
    return np.sum(k * v for k, v in energy_hist.items()) / obj_result.repetitions
    pass


pass

h, jr, jc = random_instance(3)
eprint('transverse fields: {}'.format(h))
eprint('row j fields: {}'.format(jr))
eprint('column j fields: {}'.format(jc))

circuit = cirq.Circuit()
circuit.append(one_step(h, jr, jc, 0.1, 0.2, 0.3))
eprint(circuit)

simulator = cirq.google.XmonSimulator()
circuit = cirq.Circuit()
circuit.append(one_step(h, jr, jc, 0.1, 0.2, 0.3))
circuit.append(cirq.measure(*qubits, key='x'))
results = simulator.run(circuit, repetitions=100, qubit_order=qubits)
eprint(results.histogram(key='x'))
eprint(results.histogram(key='x', fold_func=energy_func(3, h, jr, jc)))

print('>>> Value of the objective function {}'.format(obj_func(results)))

circuit = cirq.Circuit()
alpha = cirq.Symbol('alpha')
beta = cirq.Symbol('beta')
gamma = cirq.Symbol('gamma')
circuit.append(one_step(h, jr, jc, alpha, beta, gamma))
circuit.append(cirq.measure(*qubits, key='x'))
resolver = cirq.ParamResolver({'alpha': 0.1, 'beta': 0.3, 'gamma': 0.7})
resolved_circuit = circuit.with_parameters_resolved_by(resolver)
print('Resulting parameterized symbolic circuit:\n{}'.format(resolved_circuit))

spinner = Halo(text='Sweeping across the circuit', spinner='simpleDots')
spinner.start()

sweep_size = 10
sweep = (cirq.Linspace(key='alpha', start=0.0, stop=1.0, length=10)
         * cirq.Linspace(key='beta', start=0.0, stop=1.0, length=10)
         * cirq.Linspace(key='gamma', start=0.0, stop=1.0, length=10))
results = simulator.run_sweep(circuit, params=sweep, repetitions=100)

obj_min = None
min_params = None

spinner.stop_and_persist('✓')
spinner.start('Solving Objective Function')

for result in results:
    value = obj_func(result)
    if obj_min is None or value < obj_min:
        obj_min = value
        min_params = result.params
        pass
    pass
pass

spinner.stop_and_persist('✓')

print('Minimum objective value is {}.'.format(obj_min))
