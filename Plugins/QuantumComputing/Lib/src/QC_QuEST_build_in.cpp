#include "../header/QC_QuEST_build_in.h"
namespace QC
{
	void QC_QuEST_build_in::initialize_register(int number_qubits)
	{
		this->env = getQuESTEnv();
		this->qbitregister = createQureg(number_qubits);
		initZeroState(this->qbitregister);
		this->register_size = number_qubits;
	}


	void QC_QuEST_build_in::prepare_Statevector(QC::QC_Statevector Statevector)
	{
		std::vector<QC::QC_Amplitude> Amps = Statevector.get_amplitudes();
		qcomp* amp = new qcomp [ Amps.size()];
		for (int j = 0; j < Amps.size(); j++) amp[j]=getQcomp(Amps[j].real_amp,Amps[j].imag_amp );
		setQuregAmps(this->qbitregister, 0, amp, 1);
	}

	void QC_QuEST_build_in::perform_circuit(QC::QGateSequence algorithm)
	{
		std::vector<QC::QGate> Sequence = algorithm.get_sequence();
		for (int i = 0; i < Sequence.size(); i++)
		{
			QC::QGate Gate = Sequence[i];
			if (Gate.get_term() == "X")
			{
				applyPauliX(qbitregister, Gate.get_qubits_acting_on()[0]);
			}
			else if (Gate.get_term() == "Y")
			{
				applyPauliY(qbitregister, Gate.get_qubits_acting_on()[0]);
			}
			else if (Gate.get_term() == "Z")
			{
				applyPauliZ(qbitregister, Gate.get_qubits_acting_on()[0]);
			}
			else if (Gate.get_term() == "RX")
			{
				applyRotateX(qbitregister, Gate.get_qubits_acting_on()[0], Gate.get_parameters()[0]);
			}
			else if (Gate.get_term() == "RY")
			{
				applyRotateY(qbitregister, Gate.get_qubits_acting_on()[0], Gate.get_parameters()[0]);
			}
			else if (Gate.get_term() == "RZ")
			{
				applyRotateZ(qbitregister, Gate.get_qubits_acting_on()[0], Gate.get_parameters()[0]);
			}
			else if (Gate.get_term() == "H")
			{
				applyHadamard(qbitregister, Gate.get_qubits_acting_on()[0]);
			}
			else if (Gate.get_term() == "CNOT")
			{
				int* target = new int(Gate.get_qubits_acting_on()[1]);
				applyControlledMultiQubitNot(qbitregister, Gate.get_qubits_acting_on()[0], target,1);
			}
			else if (Gate.get_term() == "CRX")
			{
				applyControlledRotateX(qbitregister, Gate.get_qubits_acting_on()[0], Gate.get_qubits_acting_on()[1], Gate.get_parameters()[0]);
			}
			else if (Gate.get_term() == "CRY")
			{
				applyControlledRotateY(qbitregister, Gate.get_qubits_acting_on()[0], Gate.get_qubits_acting_on()[1], Gate.get_parameters()[0]);
			}
			else if (Gate.get_term() == "CRZ")
			{
				applyControlledRotateZ(qbitregister, Gate.get_qubits_acting_on()[0], Gate.get_qubits_acting_on()[1], Gate.get_parameters()[0]);
			}
			else if (Gate.get_term() == "CCNOT")
			{
				CompMatr2 u;
				u.real[0][0] = 0.0; u.real[0][1] = 1.0; u.real[1][0] = 1.0; u.real[1][1] = 0.0;
				u.imag[0][0] = 0.0; u.imag[0][1] = 0.0; u.imag[1][0] = 0.0; u.imag[1][1] = 0.0;
				int controll_qbits[2] = { Gate.get_qubits_acting_on()[0], Gate.get_qubits_acting_on()[1] };
				multiControlledUnitary(qbitregister, controll_qbits, 2, Gate.get_qubits_acting_on()[2], u);
			}
			else if (Gate.get_term() == "CP")
			{
				controlledPhaseShift(qbitregister, Gate.get_qubits_acting_on()[0], Gate.get_qubits_acting_on()[1], Gate.get_parameters()[0]);
			}
			else if (Gate.get_term() == "SWAP")
			{
				swapGate(qbitregister, Gate.get_qubits_acting_on()[0], Gate.get_qubits_acting_on()[1]);
			}
		}
	}

	long long int QC_QuEST_build_in::measure_with_collaps(std::vector<int> qubits_to_measure)
	{
		std::vector<string> measure_result;
		int single_qbit_result = applyMultiQubitMeasurement(qbitregister, qubits_to_measure,qubits_to_measure.size());
		measure_result = std::to_string(single_qbit_result) + measure_result;
	

		long long int measure_bin = std::stoll(measure_result);
		return measure_bin;
	}

	long long int QC_QuEST_build_in::measure_without_collaps(std::vector<int> qubits_to_measure)
	{
		Qureg copy_Register = createCloneQureg(this->qbitregister);
		std::string measure_result;
		for (int qubit_index = 0; qubit_index < qubits_to_measure.size(); qubit_index++)
		{
			int single_qbit_result = applyMultiQubitMeasurement(copy_Register, qubits_to_measure[qubit_index]);
			measure_result = std::to_string(single_qbit_result) + measure_result;
		}
		destroyQureg(copy_Register);

		long long int measure_bin = std::stoll(measure_result);
		return measure_bin;
	}

	QC::QC_Statevector QC_QuEST_build_in::get_partial_statevector(std::vector<int> amps)
	{
		QC::QC_Statevector dumpedvec;
		dumpedvec.set_number_qubits(amps[amps.size() - 1]);

		for (int i = 0; i < amps.size(); i++)
		{
			qcomp amp = getAmp(qbitregister, amps[i]);
			dumpedvec.add_amplitude(amps[i], amp.real, amp.imag);
		}

		return dumpedvec;
	}

	double QC_QuEST_build_in::get_one_qubit_probabilitie(int qubit)
	{
		qreal prop = calcProbOfOutcome(qbitregister, qubit, 1);
		return prop;
	}

	double QC_QuEST_build_in::get_pauli_expectation_value(QC::Pauli_String pauli_string)
	{
		std::vector<int> qubits_acting_on = pauli_string.get_corresponding_Qubits();
		std::vector<std::string> pauli_as_double = pauli_string.get_PauliOperations();
		std::vector<pauliOpType> paulis;
		for (int i = 0; i < pauli_as_double.size(); i++)
		{
			if (pauli_as_double[i] == "I") { paulis.push_back(PAULI_I); }
			else if (pauli_as_double[i] == "X") { paulis.push_back(PAULI_X); }
			else if (pauli_as_double[i] == "Y") { paulis.push_back(PAULI_Y); }
			else if (pauli_as_double[i] == "Z") { paulis.push_back(PAULI_Z); }
		}

		Qureg tempregister = createQureg(this->register_size, env);
		qreal prop = calcExpecPauliProd(qbitregister, &qubits_acting_on[0], &paulis[0], qubits_acting_on.size(), tempregister);
		destroyQureg(tempregister, env);
		return prop;
	}

	double QC_QuEST_build_in::get_pauli_expectation_value_with_measurement(QC::Pauli_String pauli_string, int measurements)
	{
		const double pi = 3.1415926535;
		Qureg copy_Register_Rotation = createCloneQureg(this->qbitregister, this->env);
		std::vector<int> qubits_acting_on = pauli_string.get_corresponding_Qubits();
		std::vector<std::string> paulis = pauli_string.get_PauliOperations();
		for (int i = 0; i < paulis.size(); i++)
		{
			if (paulis[i] == "X") { rotateY(copy_Register_Rotation, qubits_acting_on[i], (-pi / 2.0)); }
			else if (paulis[i] == "Y") { rotateX(copy_Register_Rotation, qubits_acting_on[i], (pi / 2.0)); }
		}

		double expectation_val_sum = 0;

		for (int i = 0; i < measurements; i++)
		{
			Qureg temp_Register = createCloneQureg(copy_Register_Rotation, this->env);
			double expectation_val = 1;
			for (int j = 0; j < paulis.size(); j++)
			{
				if (paulis[j] != "I")
				{
					int val = measure(temp_Register, qubits_acting_on[j]);
					if (val == 0) 
					{
						expectation_val *= 1;
					}
					else
					{
						expectation_val *= -1;
					}
				}
			}
			expectation_val_sum += expectation_val;
			destroyQureg(temp_Register, env);
		}
		destroyQureg(copy_Register_Rotation, env);
		return expectation_val_sum / (double)measurements;
	}


	void QC_QuEST_build_in::destroy_register()
	{
		destroyQureg(qbitregister);
		destroyQuESTEnv(env);
	}
}