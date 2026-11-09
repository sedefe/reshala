#pragma once

#include <memory>

#include "reshala/io/utils.h"
#include "reshala/model/milp_model.h"

namespace reshala {

class Io;
class Presolver;
class DualSimplex;
class MipTracker;
class HeuristicManager;
class Cutter;
class BnbSolver;
class MilpSolver;

class Reshala {
   public:
    Reshala();
    ~Reshala();

    FileReadStatus Read(const char* path);
    const Io& GetIo() const { return *io_; }

    const MilpModel& GetModel() const { return model_; }
    MilpModel& GetModel() { return model_; }

    const Presolver& GetPresolver() const { return *presolver_; }
    Presolver& GetPresolver() { return *presolver_; }

    const DualSimplex& GetDs() const { return *ds_; }
    DualSimplex& GetDs() { return *ds_; }

    const MipTracker& GetMip() const { return *mip_tracker_; }
    MipTracker& GetMip() { return *mip_tracker_; }

    const HeuristicManager& GetHeurMng() const { return *heur_manager_; }
    HeuristicManager& GetHeurMng() { return *heur_manager_; }

    const Cutter& GetCutter() const { return *cutter_; }
    Cutter& GetCutter() { return *cutter_; }

    const BnbSolver& GetBnb() const { return *bnb_; }
    BnbSolver& GetBnb() { return *bnb_; }

    Solution Solve();

    void PrintModel(std::ostream& os, const MilpModel& model) const;
    void PrintSolution(std::ostream& os, const std::vector<Scalar>& x) const;
    void PrintStats(std::ostream& os) const;

   private:
    MilpModel model_;

    std::unique_ptr<Io> io_;
    std::unique_ptr<Presolver> presolver_;
    std::unique_ptr<DualSimplex> ds_;
    std::unique_ptr<MipTracker> mip_tracker_;
    std::unique_ptr<HeuristicManager> heur_manager_;
    std::unique_ptr<Cutter> cutter_;
    std::unique_ptr<BnbSolver> bnb_;
    std::unique_ptr<MilpSolver> milp_;
};

}  // namespace reshala
