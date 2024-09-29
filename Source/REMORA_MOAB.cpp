
#include <AMReX_Utility.H>
#include "REMORA.H"

using namespace amrex;

void
REMORA::InitMOABMesh()
{
    moab::Interface * mbi = new( std::nothrow ) moab::Core;
    moab::ParallelComm* pco = new moab::ParallelComm( mbi,  ParallelContext::CommunicatorSub() );
    // Number of cells in this "domain" at this level
    std::vector<int> n_cells;

    int lev = 0;
    int which_subdomain = 0;
    int nblocks = grids[lev].size();

    // We only do single-level mesh instances
    int flev = 1; //max_level;

    Box subdomain;
    if (lev == 0) {
        subdomain = geom[lev].Domain();
    } else {
        subdomain = boxes_at_level[lev][which_subdomain];
    }

    int nx = subdomain.length(0);
    int ny = subdomain.length(1);
    int nz = subdomain.length(2);

    n_cells.push_back(nx);
    n_cells.push_back(ny);
    n_cells.push_back(nz);

    int num_vertices   = (nx+1)*(ny+1)*(nz+1);
    int num_cells = nx * ny * nz;

    std::vector<double> coords;// just the center of cells?
    long unsigned goffset = 0;
    long unsigned glen    = 0;
    int icell = 0;
    for (int i = 0; i < grids[lev].size(); ++i) {
        auto box = grids[lev][i];
        if (subdomain.contains(box)) {
            RealBox gridloc = RealBox(grids[lev][i], geom[lev].CellSize(), geom[lev].ProbLo());

            for (auto k1 = 0; k1 < grids[lev][i].length(0); ++k1) {
                for (auto k2 = 0; k2 < grids[lev][i].length(1); ++k2) {
                    for (auto k3 = 0; k3 < grids[lev][i].length(2); ++k3) {
                        coords.push_back(gridloc.lo(0)+geom[lev].CellSize(0)*static_cast<Real>(k1));
                        coords.push_back(gridloc.lo(1)+geom[lev].CellSize(1)*static_cast<Real>(k2));
                        coords.push_back(gridloc.lo(2)+geom[lev].CellSize(2)*static_cast<Real>(k3));
                        icell++;
                    }
                }
            }
            goffset += glen;
            glen = grids[lev][i].length(0)*grids[lev][i].length(1)*grids[lev][i].length(2);
        }
    }
    moab::Range cell_centers;
    moab::ErrorCode rval = mbi->create_vertices( &coords[0], icell, cell_centers );
    moab::EntityHandle moabSet;
    mbi->create_meshset( moab::MESHSET_SET, moabSet );
    mbi->add_entities(moabSet, cell_centers);

    std::string pw( "PARALLEL=WRITE_PART" );
    rval = mbi->write_file("cell_centers.h5m", 0, pw.c_str(), &moabSet, 1);


    return;
}
