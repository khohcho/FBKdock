# Third-Party Notices

FBKdock's own source code is licensed under the GNU General Public License
v3.0 or later — see [LICENSE](LICENSE). The components listed below are
bundled in this repository for convenience and remain under their original
licenses. They are aggregated with, not covered by, FBKdock's GPLv3 license.

| Component | Location | License | Upstream |
|---|---|---|---|
| Vina-GPU 2.1 (binaries + OpenCL kernel sources) | `Vina-GPU/` | Apache License 2.0 | <https://github.com/DeltaGroupNJUPT/Vina-GPU-2.1> |
| AutoDock Vina (vina, vina_split) | `AutodockVina/AutodockVina/` | Apache License 2.0 | <https://github.com/ccsb-scripps/AutoDock-Vina> |
| Boost C++ Libraries (`boost_*.dll`) | `Vina-GPU/*/`, `AutodockVina/AutodockVina/` | Boost Software License 1.0 | <https://www.boost.org/users/license.html> |
| OpenCL (Khronos headers, `OpenCL.dll`) | `Vina-GPU/*/OpenCL/` | MIT (headers); the DLL is subject to its vendor's terms | <https://www.khronos.org/opencl/> |

ADFRsuite is **not** bundled with FBKdock. It is an external dependency that
each user downloads separately from <https://ccsb.scripps.edu/adfr/downloads/>
under its own license (<https://ccsb.scripps.edu/adfr/license/>).
