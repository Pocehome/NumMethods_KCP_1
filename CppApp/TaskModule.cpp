#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "TestTask.h"
#include "KSR_Task.h"

namespace py = pybind11;

PYBIND11_MODULE(KSR1_Module, m) {
    py::class_<TestTask::FinalReference>(m, "FinalReference")
        .def_readwrite("ITERATIONS_COUNT", &TestTask::FinalReference::ITERATIONS_COUNT)
        .def_readwrite("DISTANCE_B_LAST_POINT", &TestTask::FinalReference::DISTANCE_B_LAST_POINT)
        .def_readwrite("MAX_ERROR", &TestTask::FinalReference::MAX_ERROR)
        .def_readwrite("MIN_ERROR", &TestTask::FinalReference::MIN_ERROR)
        .def_readwrite("STEP_DOUBLING_COUNT", &TestTask::FinalReference::STEP_DOUBLING_COUNT)
        .def_readwrite("STEP_REDUCTION_COUNT", &TestTask::FinalReference::STEP_REDUCTION_COUNT)
        .def_readwrite("MAX_STEP_WITH_X", &TestTask::FinalReference::MAX_STEP_WITH_X)
        .def_readwrite("MIN_STEP_WITH_X", &TestTask::FinalReference::MIN_STEP_WITH_X)
        .def_readwrite("MAX_DISTANCE_U_V", &TestTask::FinalReference::MAX_DISTANCE_U_V)
        .def_readwrite("LAST_X", &TestTask::FinalReference::LAST_X)
        .def_readwrite("LAST_V", &TestTask::FinalReference::LAST_V)
        .def_readwrite("IS_INF", &TestTask::FinalReference::IS_INF);

    //py::class_<TestTask>(m, "TestTask")
    //    .def(py::init<double, double, double, double, double, int, double>())
    //    .def("Solve_Without_Error_Control", &TestTask::Solve_Without_Error_Control)
    //    .def("Solve_With_Error_Control", &TestTask::Solve_With_Error_Control)
    //    .def("get_table_information", &TestTask::get_table_information)
    //    .def("get_final_reference", &TestTask::get_reference);

    py::class_<KSR_Task>(m, "KSR_Task")
        .def(py::init<double, double, double, double, double, int, double>())
        .def("Solve_Without_Error_Control", &KSR_Task::Solve_Without_Error_Control)
        .def("Solve_With_Error_Control", &KSR_Task::Solve_With_Error_Control)
        .def("get_table_information", &KSR_Task::get_table_information)
        .def("get_final_reference", &KSR_Task::get_reference)
        .def("set_params", &KSR_Task::set_params)
        .def("get_ERRORS_LIST", &KSR_Task::get_ERRORS_LIST);
}
