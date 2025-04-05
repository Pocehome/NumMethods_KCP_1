#pragma once
#include "TestTask.h"
#include <vector>


class SecondTask:public TestTask
{

protected:

    std::vector<double> f(const double& X, const std::vector<double>& V);
    std::vector<double> find_K1(const double& X, const std::vector<double>& V);													                // �������, ������� ������� K1 �� ������
    std::vector<double> find_K2(const double& X, const std::vector<double>& V, const double& STEP, const std::vector<double>& K1);				// ...K2
    std::vector<double> find_K3(const double& X, const std::vector<double>& V, const double& STEP, const std::vector<double>& K2);				// ...K3
    std::vector<double> find_K4(const double& X, const std::vector<double>& V, const double& STEP, const std::vector<double>& K3);              // ...K4
    void make_Step(double& X, std::vector<double>& V, const double& STEP);

    int control_Error(double& X, std::vector<double>& V, double& X_EXTRA, std::vector<double>& V_EXTRA, double& OLD_X,
                       std::vector<double>& OLD_V, double& S, double& CURRENT_DOUBLING, double& CURRENT_REDUCTION);
    double alpha{}; //�������� � ��� ���������
public:

    void set_alpha(double _alpha); //�������� � ��� ���������

    SecondTask(const double& _A, const double& _B, const double& _STEP, const double& _E_ERROR,
        const double& _E_BORDER, const int& _MAX_STEPS, const double& _START_POINT)
        : TestTask(_A, _B, _STEP, _E_ERROR, _E_BORDER, _MAX_STEPS, _START_POINT) {};

    void Solve_Without_Error_Control() override;
    void Solve_With_Error_Control() override;

    // ������� ��� �������� � main , � �������� ��������, ��� ����� �� �����
    void PrintTable() override;     //����� �������� ������� �� �������, ���� ����� ������ ������ �������, ����� �� �����������
    void Write_To_File() override;  //������ X,Vi,Ui � ����
    void PrintReference() override; //����� �������� �������
};

