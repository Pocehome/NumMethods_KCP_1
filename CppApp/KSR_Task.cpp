#include "KSR_Task.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>

bool is_infin(const std::vector<double>& V, const std::vector<double>& V_EX) {
	bool flag{ false };
	if (std::isinf(V[0]) || std::isnan(V[0])) flag = true;
	if (std::isinf(V[1]) || std::isnan(V[1])) flag = true;
	if (std::isinf(V_EX[0]) || std::isnan(V_EX[0])) flag = true;
	if (std::isinf(V_EX[1]) || std::isnan(V_EX[1])) flag = true;

	return flag;
}

std::vector<double> KSR_Task::f(const double& X, const std::vector<double>& V) {
	double f_1 = V[1];
	double f_2 = -c/m*V[1] - k/m*V[0] - k_with_star/m*pow(V[0], 3);

	std::vector<double> function(2);
	function[0] = f_1;
	function[1] = f_2;

	return function;
}

std::vector<double> KSR_Task::find_K1(const double& X, const std::vector<double>& V) {
	return f(X, V);
}

std::vector<double> KSR_Task::find_K2(const double& X, const std::vector<double>& V, const double& STEP, const std::vector<double>& K1) {
	std::vector<double> new_V{ V[0] + (STEP / 2.0) * K1[0], V[1] + (STEP / 2.0) * K1[1] };
	return f(X + STEP / 2.0, new_V);
}

std::vector<double> KSR_Task::find_K3(const double& X, const std::vector<double>& V, const double& STEP, const std::vector<double>& K2) {
	std::vector<double> new_V{ V[0] + (STEP / 2.0) * K2[0], V[1] + (STEP / 2.0) * K2[1] };
	return f(X + STEP / 2.0, new_V);
}

std::vector<double> KSR_Task::find_K4(const double& X, const std::vector<double>& V, const double& STEP, const std::vector<double>& K3) {
	std::vector<double> new_V{ V[0] + STEP * K3[0], V[1] + STEP * K3[1] };
	return f(X + STEP, new_V);
}

void KSR_Task::make_Step(double& X, std::vector<double>& V, const double& STEP) {
	std::vector<double> K1 = find_K1(X, V);
	std::vector<double> K2 = find_K2(X, V, STEP, K1);
	std::vector<double> K3 = find_K3(X, V, STEP, K2);
	std::vector<double> K4 = find_K4(X, V, STEP, K3);

	X += STEP;
	V[0] += (STEP / 6.0) * (K1[0] + 2.0 * K2[0] + 2.0 * K3[0] + K4[0]);
	V[1] += (STEP / 6.0) * (K1[1] + 2.0 * K2[1] + 2.0 * K3[1] + K4[1]);
}

int KSR_Task::control_Error(double& X, std::vector<double>& V, double& X_EXTRA, std::vector<double>& V_EXTRA, double& OLD_X, std::vector<double>& OLD_V, double& S, double& CURRENT_DOUBLING, double& CURRENT_REDUCTION) {
	if (abs(S) < parametrs.E_ERROR / pow(2.0, 5)) { // Условие для удвоения
		parametrs.STEP *= 2.0;
		++reference.STEP_DOUBLING_COUNT;
		++CURRENT_DOUBLING;

		return 0;
	}
	else if (abs(S) > parametrs.E_ERROR) {          //Если условие для деления шага выполнилось, то возвращаемся назад,
		bool FLAG_TO_EXIT = false;
		double h = parametrs.STEP;  				// выполняем расчёт с половинным шагом, если для него погрешность снова 
		while (!FLAG_TO_EXIT) {
			// оказалась большой, то снова возвращаемся и делим шаг, и так пока погрешность не будет допустимой
			h /= 2.0;
			++CURRENT_REDUCTION;
			++reference.STEP_REDUCTION_COUNT;

			V = OLD_V;
			X = OLD_X;
			V_EXTRA = V;
			X_EXTRA = X;

			for (int i = 0; i < 2; ++i) {
				make_Step(X_EXTRA, V_EXTRA, h / 2.0);
			}
			make_Step(X, V, h);

			double s_1 = (V_EXTRA[0] - V[0]) / (pow(2.0, 4.0) - 1.0);
			double s_2 = (V_EXTRA[1] - V[1]) / (pow(2.0, 4.0) - 1.0);
			S = sqrt(s_1 * s_1 + s_2 * s_2);
			if (abs(S) <= parametrs.E_ERROR) {
				FLAG_TO_EXIT = true;
				parametrs.STEP = h;
			}
		}

		return 1;
	}
}

void KSR_Task::set_params(double _k, double _k_with_star, double _c, double _m) {
	k = _k;
	k_with_star = _k_with_star;
	c = _c;
	m = _m;
}

void KSR_Task::Solve_Without_Error_Control() {
	double X = parametrs.A;
	std::vector<double> V = { parametrs.START_POINT_FOR_U, 0.0 }; //Далее V[0] это U - искомая функция, V[1] это U' - её производная 

	std::vector<double> TABLE_ROW1 = { 0.0, X, V[0],V[1], parametrs.STEP }; // Здесь и далее - строка итоговой таблицы в виде (5 элементов)
	TABLE_INFORMATION.emplace_back(TABLE_ROW1);
	// { i; X_i; V0_i; V1_i; STEP_i }

	double OLD_X = X; //В переменых OLD храним значения с последного шага
	std::vector<double> OLD_V = V;

	bool FLAG_TO_EXIT = false;
	for (int i = 1; i <= parametrs.MAX_STEPS; ++i)
	{
		make_Step(X, V, parametrs.STEP);

		std::vector<double> tmp{ 1.0,1.0 };
		if (is_infin(V, tmp)) {
			reference.IS_INF = true;
			break;
		}


		if (X > parametrs.B) {						// Если X вышел за правую границу, возвращаемся на шаг назад и делаем шаг,                                        
			X = OLD_X;								//  равный разнице правой границы и последней точкой X и заканчиваем интегрирование
			V = OLD_V;

			parametrs.STEP = parametrs.B - OLD_X;
			make_Step(X, V, parametrs.STEP);

			FLAG_TO_EXIT = true;                       //Если X совпал с правой границей, заканчиваем интегрирование 
		}

		std::vector<double> TABLE_ROW = { static_cast<double>(i), X, V[0],V[1], parametrs.STEP };
		TABLE_INFORMATION.emplace_back(TABLE_ROW);
		++reference.ITERATIONS_COUNT;

		if (FLAG_TO_EXIT)
			break;
		else {
			OLD_X = X;
			OLD_V = V;
		}


	}
}

void KSR_Task::Solve_With_Error_Control() {
	double X = parametrs.A;
	std::vector<double> V = { parametrs.START_POINT_FOR_U,0.0 };

	std::vector<double> TABLE_ROW1 = { 0.0, X, V[0], V[1], V[0], V[1], 0.0, 0.0, 0.0, parametrs.STEP, 0.0, 0.0 }; // Здесь и далее - строка итоговой таблицы в виде (12 элементов)
	TABLE_INFORMATION.emplace_back(TABLE_ROW1);									// { i; X_i; V[0]_i; V[1]_i; V[0]_i^; V[1]_i^; V[0]_i-V[0]_i^; V[1]_i-V[1]_i^; ОЛП(S); STEP_i; Кол-во делений; Кол-во удвоений  }

	double OLD_X = X;
	std::vector<double> OLD_V = V;

	for (int i = 1; i <= parametrs.MAX_STEPS; ++i)
	{
		double CURRENT_DOUBLING{};				// Кол-во удвоений на текущем шаге
		double CURRENT_REDUCTION{};				// Кол-во делений на текущем шаге

		std::vector<double> V_EXTRA = V;         // V^ для половинного шага
		double X_EXTRA = X;						 // X^ для половинного шага

		for (int j = 0; j < 2; ++j) {
			make_Step(X_EXTRA, V_EXTRA, parametrs.STEP / 2.0);
		}
		make_Step(X, V, parametrs.STEP);

		if (is_infin(V, V_EXTRA)) {
			reference.IS_INF = true;
			break;
		}

		double s_1 = (V_EXTRA[0] - V[0]) / (pow(2.0, 4.0) - 1.0);
		double s_2 = (V_EXTRA[1] - V[1]) / (pow(2.0, 4.0) - 1.0);
		double S = sqrt(s_1 * s_1 + s_2 * s_2);

		double H = parametrs.STEP;
		int step_control = control_Error(X, V, X_EXTRA, V_EXTRA, OLD_X, OLD_V, S, CURRENT_DOUBLING, CURRENT_REDUCTION); // Непосредственно сам контроль ЛП, подробнее см. в реализации функции

		double step;
		if (step_control == 0) {
			step = H;
		}
		else {
			step = parametrs.STEP;
		}
		double LP = abs(S * pow(2.0, 4));
		ERRORS_LIST.emplace_back(LP);

		if (X < parametrs.B - parametrs.E_BORDER) {
			std::vector<double> TABLE_ROW = { static_cast<double>(i), X, V[0], V[1], V_EXTRA[0], V_EXTRA[1], V[0]-V_EXTRA[0], V[1]-V_EXTRA[1], LP, step, CURRENT_REDUCTION, CURRENT_DOUBLING };
			TABLE_INFORMATION.emplace_back(TABLE_ROW);

			++reference.ITERATIONS_COUNT;
			STEPS_and_Xs.emplace_back(std::make_pair(parametrs.STEP, X));
		}

		bool EXIT_FROM_FOR = false;
		if (x_in_border(parametrs.B, X, parametrs.E_BORDER)) {                       //Проверка на попадание в окрестность правой границы по X. 
			std::vector<double> TABLE_ROW = { static_cast<double>(i), X, V[0], V[1], V_EXTRA[0], V_EXTRA[1], V[0]-V_EXTRA[0], V[1]-V_EXTRA[1], LP, step, CURRENT_REDUCTION, CURRENT_DOUBLING };
			TABLE_INFORMATION.emplace_back(TABLE_ROW);

			++reference.ITERATIONS_COUNT;
			STEPS_and_Xs.emplace_back(std::make_pair(parametrs.STEP, X));

			EXIT_FROM_FOR = true;												   //Если X попал в окрестность, завершаем интегрирование, выходя из for по флагу
		}

		if (X > parametrs.B) {						// Если X вышел за правую границу, возвращаемся на шаг назад и делаем шаг,                                        
			X = OLD_X;								//  равный разнице правой границы и последней точкой X и заканчиваем интегрирование
			V = OLD_V;

			parametrs.STEP = parametrs.B - OLD_X;
			make_Step(X, V, parametrs.STEP);

			std::vector<double> TABLE_ROW = { static_cast<double>(i), X, V[0],V[1], V_EXTRA[0],V_EXTRA[1], V[0] - V_EXTRA[0],V[1] - V_EXTRA[1], LP, parametrs.STEP, CURRENT_REDUCTION, CURRENT_DOUBLING };
			TABLE_INFORMATION.emplace_back(TABLE_ROW);
			++reference.ITERATIONS_COUNT;

			EXIT_FROM_FOR = true;                        //Если X совпал с правой границей, заканчиваем интегрирование 
		}

		if (EXIT_FROM_FOR) {
			break;
		}
		else {
			OLD_X = X;
			OLD_V = V;
		}
	}

	reference.DISTANCE_B_LAST_POINT = parametrs.B - X;
	find_Max_Step();
	find_Min_Step();
	find_Max_Error();
	find_Min_Error();
	find_Last_X();
	find_Last_V();
}

std::list<double> KSR_Task::get_ERRORS_LIST() {
	return ERRORS_LIST;
}

void KSR_Task::PrintTable() {
	if (TABLE_INFORMATION.front().size() == 5) {
		std::cout << std::setw(30) << "i" << std::setw(30) << "X" << std::setw(30) << "V[0]" << std::setw(30) << "V[1]" << std::setw(30) << "Step" << std::endl;
		for (auto it_list = TABLE_INFORMATION.begin(); it_list != TABLE_INFORMATION.end(); ++it_list) {
			for (size_t i = 0; i < 5; ++i) {
				if (i == 1) {
					std::cout << std::setw(30) << (*it_list)[i];
				}
				else {
					std::cout << std::setprecision(16) << std::setw(30) << (*it_list)[i];
				}
			}
			std::cout << std::endl;
		}
	}
	else {
		std::cout << std::setw(30) << "i" << std::setw(30) << "X" << std::setw(30) << "V[0]" << std::setw(30) << "V[1]" << std::setw(30) << "V[0]^" << std::setw(30) << "V[1]^"
			<< std::setw(30) << "V[0]-V[0]^" << std::setw(30) << "V[1]-V[1]^" << std::setw(30) << "S"
			<< std::setw(30) << "Step" << std::setw(30) << "Уменьшения" << std::setw(30) << "Удвоения" << std::setw(30) << std::endl;
		for (auto it_list = TABLE_INFORMATION.begin(); it_list != TABLE_INFORMATION.end(); ++it_list) {
			for (size_t i = 0; i < 12; ++i) {
				if (i == 1) {
					std::cout << std::setw(30) << (*it_list)[i];
				}
				else if (i == 10 || i == 11) {
					std::cout << std::setw(30) << static_cast<int>((*it_list)[i]);
				}
				else {
					std::cout << std::setprecision(16) << std::setw(30) << (*it_list)[i];
				}
			}
			std::cout << std::endl;
		}
	}
}

void KSR_Task::Write_To_File() {
	// { i; X_i; V[0]_i; V[1]_i; V[0]_i^; V[1]_i^; V[0]_i-V[0]_i^; V[1]_i-V[1]_i^; ОЛП(S); STEP_i; Кол-во делений; Кол-во удвоений  }
	std::ofstream outFile("points.txt");
	if (!outFile) {
		std::cerr << "Ошибка открытия файла!" << std::endl;
	}
	else {

		for (auto it_list = TABLE_INFORMATION.begin(); it_list != TABLE_INFORMATION.end(); ++it_list) {
			outFile << std::fixed << std::setprecision(20) << (*it_list)[1] << "\t" << std::setw(30)
				<< std::setprecision(20) << (*it_list)[2] << "\t" << std::setw(30)
				<< std::setprecision(20) << (*it_list)[3] << std::endl;  // Запись в формате "x v[0] v[1]"
		}
		outFile.close();
	}
}

void KSR_Task::PrintReference() {
	std::cout << std::endl;
	std::cout << "DISTANCE_B_LAST_POINT :" << reference.DISTANCE_B_LAST_POINT << std::endl;
	std::cout << "ITERATIONS_COUNT :" << reference.ITERATIONS_COUNT << std::endl;
	std::cout << "MAX_ERROR :" << reference.MAX_ERROR << std::endl;
	std::cout << "MAX_STEP_WITH_X :" << "STEP = " << reference.MAX_STEP_WITH_X.first << "  X = " << reference.MAX_STEP_WITH_X.second << std::endl;
	std::cout << "MIN_STEP_WITH_X :" << "STEP = " << reference.MIN_STEP_WITH_X.first << "  X = " << reference.MIN_STEP_WITH_X.second << std::endl;
	std::cout << "STEP_DOUBLING_COUNT :" << reference.STEP_DOUBLING_COUNT << std::endl;
	std::cout << "STEP_REDUCTION_COUNT :" << reference.STEP_REDUCTION_COUNT << std::endl;
}
