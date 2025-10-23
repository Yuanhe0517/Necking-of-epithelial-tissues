/*!
 * \file simulation.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \modified by Yuan He, yuanhe@westlake.edu.cn
 * \date 30-Nov-2023 <Modification Date 22-Oct-2025>
 * \brief Simulation class 
 */ 

#include "simulation.hpp"

namespace VMTutorial
{
  bool Simulation::run(int steps, bool topological_change, bool old_style)
  {
    //_integ.set_dt(0.12);
    //double progress = 0.0;
    for (int i = sim_step; i < sim_step + steps; i++)
    {
     
      if (topological_change)
      {
          t += _topology.T1();
      }
      //cout<<"apply"<<endl;
      _integ.apply();
      //cout<<"finish apply"<<endl;
      _sys.set_topology_change(false);

      // FIRE algorithm: not use for overdamped systems

      // if (_integ._power >0)
      // {
      //     _integ._dt = std::min(_integ._dt*1.1,_integ.dt_max);
      //     _integ._alpha *=0.99;
      // }
      // else{
      //     _integ._dt *=0.5;
      //     _integ._alpha = 0.1;

      // }


      if (_integ._meq)
      {
        //_dump.dump_junctions("1");
        //_dump.dump_cells("2");
        _integ.set_eq(false);
        //_integ._dt = 0.08;
        _integ.set_dt(0.12);
        appendValueToCSV(t);
        t=0;
        cout<<"Equil!"<<_integ._dt<<endl;
        return true;
      }

      
    }

    _integ._dt = _integ._dt*0.9;
    _integ.set_dt(_integ._dt);
    cout<<"Adapt!"<<_integ._dt<<endl;
    return false;

  }
  
    // if (this->print_freq > 0 && !old_style)
    //   progress_bar(progress, " ");
    
    //sim_step += steps;
    // if (this->print_freq > 0 && !old_style)
    //   std::cout << " --> Completed " << sim_step << " simulation steps " << std::endl;  
    
  

  void Simulation::progress_bar(double progress, const string& end_of_line)
  {
    std::cout << "[";
    int pos = static_cast<int>(round(bar_width * progress));
    for (int i = 0; i < bar_width; ++i) 
    {
      if (i < pos) 
        std::cout << "=";
      else if (i == pos) 
        std::cout << ">";
      else 
        std::cout << " ";
    }
    std::cout << "] "  << static_cast<int>(round(progress * 100.0)) << "%" << end_of_line;
    std::cout.flush();
  }

  // 函数用于将输入的数值k追加写入到已有的Excel文件中

  void Simulation::appendValueToCSV(int k) 
  {
    std::ofstream file;
    file.open("data.csv", std::ios::app);  // 以追加模式打开文件
    if (file.is_open()) {
        file << k << "\n";  // 将数值k写入文件，并换行
        file.close();
        std::cout << "数值 " << k << " 已成功追加到CSV文件中。" << std::endl;
    } else {
        std::cout << "无法打开CSV文件。" << std::endl;
    }
  }

  void export_Simulation(py::module& m)
  {
    py::class_<Simulation>(m, "Simulation")
          .def(py::init<System&,Integrate&,ForceCompute&,Topology&>())
          .def_readwrite("print_freq", &Simulation::print_freq)
          .def_readwrite("bar_width", &Simulation::bar_width)
          .def("run", &Simulation::run, py::arg("steps"), py::arg("topological_change") = true, py::arg("old_style") = false)
          .def("print_version", &Simulation::print_version);
  }  


}

PYBIND11_MODULE(vm, m)
{
  VMTutorial::export_Vec(m);
  VMTutorial::export_Box(m);
  VMTutorial::export_VertexProperty(m);
  VMTutorial::export_HEProperty(m);
  VMTutorial::export_EdgeProperty(m);
  VMTutorial::export_FaceProperty(m);
  VMTutorial::export_Vertex(m);
  VMTutorial::export_VertexCirculator(m);
  VMTutorial::export_Edge(m);
  VMTutorial::export_HalfEdge(m);
  VMTutorial::export_Face(m);
  VMTutorial::export_FaceCirculator(m);
  VMTutorial::export_Mesh(m);
  VMTutorial::export_System(m);
  VMTutorial::export_ForceCompute(m);
  VMTutorial::export_Integrate(m);
  VMTutorial::export_Topology(m);
  VMTutorial::export_Dump(m);
  VMTutorial::export_Simulation(m);
}

