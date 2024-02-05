#include "Driver.hpp"

#include <iostream>

Driver::Driver(const InputFile& input, const std::string& problem_name)
    : mesh{std::make_shared<Mesh>(input)},
      diffusion{input, mesh},
      writer{problem_name, mesh},
      _problem_name(problem_name) {
    std::cout << "+++++++++++++++++++++" << '\n';
    std::cout << "  Running deqn (Refactored by Xanonymous)  " << '\n';

#ifdef DEBUG
    std::cout << "- input file: " << problem_name << '\n';
#endif

    dt_max = input.getDouble("dt_max", 0.2);
    dt = input.getDouble("initial_dt", 0.2);

    t_start = input.getDouble("start_time", 0.0);
    t_end = input.getDouble("end_time", 10.0);

    vis_frequency = input.getInt("vis_frequency", -1);
    summary_frequency = input.getInt("summary_frequency", 1);

#ifdef DEBUG
    std::cout << "- dt_max: " << dt_max << '\n';
    std::cout << "- initial_dt: " << dt << '\n';
    std::cout << "- start_time: " << t_start << '\n';
    std::cout << "- end_time: " << t_end << '\n';
    std::cout << "- vis_frequency: " << vis_frequency << '\n';
    std::cout << "- summary_frequency: " << summary_frequency << '\n';
#endif

    std::cout << "+++++++++++++++++++++" << '\n';
    std::cout << '\n';

    /* Initial mesh dump */
    if (vis_frequency != -1) {
        writer.writeVtk(0, 0.0);
    }
}

void Driver::run() const {
    int step;
    double t_current = t_start;

    std::cout
            << "+++++++++++++++++++++" << '\n'
            << "     Run Started     " << '\n'
            << "+++++++++++++++++++++" << '\n';

    const int& total_runtime = static_cast<int>((t_end - t_start) / dt);

    if (vis_frequency == 1 && summary_frequency == 1) {
        for (step = 1; step <= total_runtime; ++step) {
            t_current += dt;

            diffusion.doCycle(dt);
            writer.writeVtk(step, t_current);

#ifdef DEBUG
            std::cout << "+ step: " << step << ", dt:   " << dt << '\n';
            const double temperature = mesh->getTotalTemperature();
            std::cout << "+\tcurrent total temperature: " << temperature << '\n';
#endif
        }

        writer.writeVisited(step - 1);
    } else {
        for (step = 1; step <= total_runtime; ++step) {
            t_current += dt;

            diffusion.doCycle(dt);

            if (step % vis_frequency == 0 && vis_frequency != -1) {
                writer.writeVtk(step, t_current);
            }

#ifdef DEBUG
            std::cout << "+ step: " << step << ", dt:   " << dt << '\n';
            if (step % summary_frequency == 0 && summary_frequency != -1) {
                const double temperature = mesh->getTotalTemperature();
                std::cout << "+\tcurrent total temperature: " << temperature << '\n';
            }
#endif
        }

        if (step % vis_frequency != 0 && vis_frequency != -1) {
            writer.writeVtk(step, t_current);
            writer.writeVisited(step);
        } else {
            writer.writeVisited(step - 1);
        }
    }

    std::cout << '\n';
    std::cout << "+++++++++++++++++++++" << '\n';
    std::cout << "   Run completete.   " << '\n';
    std::cout << "+++++++++++++++++++++" << '\n';
}
