# Build using:
# git clean -dfx && FC=lfortran cmake -DCMAKE_Fortran_FLAGS=--verbose -DLFORTRAN_RUNTIME_LIBRARY_PATH=/Users/ondrej/repos/lfortran/lfortran/src/runtime/liblfortran_runtime.dylib . && make -j
cmake_minimum_required(VERSION 3.10)

project(LFortranScipyBuild)
enable_language(Fortran)

set(CMAKE_SHARED_LIBRARY_PREFIX "") # Remove the "lib" prefix from the shared library
set(SPECFUN_DIRECTORY ${CMAKE_BINARY_DIR}/../scipy/special/)
set(MINPACK_DIRECTORY ${CMAKE_BINARY_DIR}/../scipy/optimize/)
set(FITPACK_DIRECTORY ${CMAKE_BINARY_DIR}/../scipy/interpolate/)
set(INTEGRATE_DIRECTORY ${CMAKE_BINARY_DIR}/../scipy/integrate/)




# Specify the source files explicitly for each component Specfun, Minpack, Fitpack
set(SOURCE_FILES_SPECFUN
    "../scipy/special/specfun/specfun.f"
    # Add more source files here if needed
)

set(SOURCE_FILES_MINPACK
    "../scipy/optimize/minpack/dogleg.f"
    "../scipy/optimize/minpack/dpmpar.f"
    "../scipy/optimize/minpack/enorm.f"
    "../scipy/optimize/minpack/fdjac1.f"
    "../scipy/optimize/minpack/fdjac2.f"
    "../scipy/optimize/minpack/hybrd.f"
    "../scipy/optimize/minpack/hybrd1.f"
    "../scipy/optimize/minpack/hybrj.f"
    "../scipy/optimize/minpack/hybrj1.f"
    "../scipy/optimize/minpack/lmder.f"
    "../scipy/optimize/minpack/lmder1.f"
    "../scipy/optimize/minpack/lmdif.f"
    "../scipy/optimize/minpack/lmdif1.f"
    "../scipy/optimize/minpack/lmpar.f"
    "../scipy/optimize/minpack/lmstr.f"
    "../scipy/optimize/minpack/lmstr1.f"
    "../scipy/optimize/minpack/qform.f"
    "../scipy/optimize/minpack/qrfac.f"
    "../scipy/optimize/minpack/qrsolv.f"
    "../scipy/optimize/minpack/r1mpyq.f"
    "../scipy/optimize/minpack/r1updt.f"
    "../scipy/optimize/minpack/rwupdt.f"
    "../scipy/optimize/minpack/chkder.f"
)

set(SOURCE_FILES_MINPACK_RTLIB
    "../scipy/optimize/minpack/chkder.f"
)

set(SOURCE_FILES_FITPACK
    "../scipy/interpolate/fitpack/bispeu.f"
    "../scipy/interpolate/fitpack/bispev.f"
    "../scipy/interpolate/fitpack/clocur.f"
    "../scipy/interpolate/fitpack/cocosp.f"
    "../scipy/interpolate/fitpack/concon.f"
    "../scipy/interpolate/fitpack/concur.f"
    "../scipy/interpolate/fitpack/cualde.f"
    "../scipy/interpolate/fitpack/curev.f"
    "../scipy/interpolate/fitpack/curfit.f"
    "../scipy/interpolate/fitpack/dblint.f"
    "../scipy/interpolate/fitpack/fourco.f"
    "../scipy/interpolate/fitpack/fpader.f"
    "../scipy/interpolate/fitpack/fpadno.f"
    "../scipy/interpolate/fitpack/fpadpo.f"
    "../scipy/interpolate/fitpack/fpback.f"
    "../scipy/interpolate/fitpack/fpbacp.f"
    "../scipy/interpolate/fitpack/fpbfout.f"
    "../scipy/interpolate/fitpack/fpbisp.f"
    "../scipy/interpolate/fitpack/fpbspl.f"
    "../scipy/interpolate/fitpack/fpchec.f"
    "../scipy/interpolate/fitpack/fpched.f"
    "../scipy/interpolate/fitpack/fpchep.f"
    "../scipy/interpolate/fitpack/fpclos.f"
    "../scipy/interpolate/fitpack/fpcoco.f"
    "../scipy/interpolate/fitpack/fpcons.f"
    "../scipy/interpolate/fitpack/fpcosp.f"
    "../scipy/interpolate/fitpack/fpcsin.f"
    "../scipy/interpolate/fitpack/fpcurf.f"
    "../scipy/interpolate/fitpack/fpcuro.f"
    "../scipy/interpolate/fitpack/fpcyt1.f"
    "../scipy/interpolate/fitpack/fpcyt2.f"
    "../scipy/interpolate/fitpack/fpdeno.f"
    "../scipy/interpolate/fitpack/fpdisc.f"
    "../scipy/interpolate/fitpack/fpfrno.f"
    "../scipy/interpolate/fitpack/fpgivs.f"
    "../scipy/interpolate/fitpack/fpgrdi.f"
    "../scipy/interpolate/fitpack/fpgrpa.f"
    "../scipy/interpolate/fitpack/fpgrre.f"
    "../scipy/interpolate/fitpack/fpgrsp.f"
    "../scipy/interpolate/fitpack/fpinst.f"
    "../scipy/interpolate/fitpack/fpintb.f"
    "../scipy/interpolate/fitpack/fpknot.f"
    "../scipy/interpolate/fitpack/fpopdi.f"
    "../scipy/interpolate/fitpack/fpopsp.f"
    "../scipy/interpolate/fitpack/fporde.f"
    "../scipy/interpolate/fitpack/fppara.f"
    "../scipy/interpolate/fitpack/fppasu.f"
    "../scipy/interpolate/fitpack/fpperi.f"
    "../scipy/interpolate/fitpack/fppocu.f"
    "../scipy/interpolate/fitpack/fppogr.f"
    "../scipy/interpolate/fitpack/fppola.f"
    "../scipy/interpolate/fitpack/fprank.f"
    "../scipy/interpolate/fitpack/fprati.f"
    "../scipy/interpolate/fitpack/fpregr.f"
    "../scipy/interpolate/fitpack/fprota.f"
    "../scipy/interpolate/fitpack/fprppo.f"
    "../scipy/interpolate/fitpack/fprpsp.f"
    "../scipy/interpolate/fitpack/fpseno.f"
    "../scipy/interpolate/fitpack/fpspgr.f"
    "../scipy/interpolate/fitpack/fpsphe.f"
    "../scipy/interpolate/fitpack/fpsuev.f"
    "../scipy/interpolate/fitpack/fpsurf.f"
    "../scipy/interpolate/fitpack/fpsysy.f"
    "../scipy/interpolate/fitpack/fptrnp.f"
    "../scipy/interpolate/fitpack/fptrpe.f"
    "../scipy/interpolate/fitpack/insert.f"
    "../scipy/interpolate/fitpack/parcur.f"
    "../scipy/interpolate/fitpack/parder.f"
    "../scipy/interpolate/fitpack/pardeu.f"
    "../scipy/interpolate/fitpack/parsur.f"
    "../scipy/interpolate/fitpack/percur.f"
    "../scipy/interpolate/fitpack/pogrid.f"
    "../scipy/interpolate/fitpack/polar.f"
    "../scipy/interpolate/fitpack/profil.f"
    "../scipy/interpolate/fitpack/regrid.f"
    "../scipy/interpolate/fitpack/spalde.f"
    "../scipy/interpolate/fitpack/spgrid.f"
    "../scipy/interpolate/fitpack/sphere.f"
    "../scipy/interpolate/fitpack/splder.f"
    "../scipy/interpolate/fitpack/splev.f"
    "../scipy/interpolate/fitpack/splint.f"
    "../scipy/interpolate/fitpack/sproot.f"
    "../scipy/interpolate/fitpack/surev.f"
    "../scipy/interpolate/fitpack/surfit.f"
    "../scipy/interpolate/fitpack/pardtc.f"
    "../scipy/interpolate/fitpack/evapol.f"
)

SET(SOURCE_FILES_FITPACK_RTLIB
    "../scipy/interpolate/fitpack/fpcuro.f"
)

set(SOURCE_FILES_AMOS
    "../scipy/special/amos/dgamln.f"
    "../scipy/special/amos/zabs.f"
    "../scipy/special/amos/zacai.f"
    "../scipy/special/amos/zacon.f"
    "../scipy/special/amos/zairy.f"
    "../scipy/special/amos/zasyi.f"
    "../scipy/special/amos/zbesh.f"
    "../scipy/special/amos/zbesi.f"
    "../scipy/special/amos/zbesj.f"
    "../scipy/special/amos/zbesk.f"
    "../scipy/special/amos/zbesy.f"
    "../scipy/special/amos/zbinu.f"
    "../scipy/special/amos/zbiry.f"
    "../scipy/special/amos/zbknu.f"
    "../scipy/special/amos/zbuni.f"
    "../scipy/special/amos/zbunk.f"
    "../scipy/special/amos/zdiv.f"
    "../scipy/special/amos/zexp.f"
    "../scipy/special/amos/zkscl.f"
    "../scipy/special/amos/zlog.f"
    "../scipy/special/amos/zmlri.f"
    "../scipy/special/amos/zmlt.f"
    "../scipy/special/amos/zrati.f"
    "../scipy/special/amos/zs1s2.f"
    "../scipy/special/amos/zseri.f"
    "../scipy/special/amos/zshch.f"
    "../scipy/special/amos/zsqrt.f"
    "../scipy/special/amos/zuchk.f"
    "../scipy/special/amos/zunhj.f"
    "../scipy/special/amos/zuni1.f"
    "../scipy/special/amos/zuni2.f"
    "../scipy/special/amos/zunik.f"
    "../scipy/special/amos/zunk1.f"
    "../scipy/special/amos/zunk2.f"
    "../scipy/special/amos/zuoik.f"
    "../scipy/special/amos/zwrsk.f"
    "../scipy/special/mach/d1mach.f"
    "../scipy/special/mach/i1mach.f"
    "../scipy/special/mach/xerror.f"
)

SET(SOURCE_FILES_AMOS_RTLIB
    "../scipy/special/amos/zairy.f"
)

set(SOURCE_FILES_MACH
    "../scipy/special/mach/d1mach.f"
    "../scipy/special/mach/i1mach.f"
    "../scipy/special/mach/xerror.f"
)

set(SOURCE_FILES_CDFLIB
    "../scipy/special/cdflib/algdiv.f"
    "../scipy/special/cdflib/alngam.f"
    "../scipy/special/cdflib/alnrel.f"
    "../scipy/special/cdflib/apser.f"
    "../scipy/special/cdflib/basym.f"
    "../scipy/special/cdflib/bcorr.f"
    "../scipy/special/cdflib/betaln.f"
    "../scipy/special/cdflib/bfrac.f"
    "../scipy/special/cdflib/bgrat.f"
    "../scipy/special/cdflib/bpser.f"
    "../scipy/special/cdflib/bratio.f"
    "../scipy/special/cdflib/brcmp1.f"
    "../scipy/special/cdflib/brcomp.f"
    "../scipy/special/cdflib/bup.f"
    "../scipy/special/cdflib/cdfbet.f"
    "../scipy/special/cdflib/cdfbin.f"
    "../scipy/special/cdflib/cdfchi.f"
    "../scipy/special/cdflib/cdfchn.f"
    "../scipy/special/cdflib/cdff.f"
    "../scipy/special/cdflib/cdffnc.f"
    "../scipy/special/cdflib/cdfgam.f"
    "../scipy/special/cdflib/cdfnbn.f"
    "../scipy/special/cdflib/cdfnor.f"
    "../scipy/special/cdflib/cdfpoi.f"
    "../scipy/special/cdflib/cdft.f"
    "../scipy/special/cdflib/cdftnc.f"
    "../scipy/special/cdflib/cumbet.f"
    "../scipy/special/cdflib/cumbin.f"
    "../scipy/special/cdflib/cumchi.f"
    "../scipy/special/cdflib/cumchn.f"
    "../scipy/special/cdflib/cumf.f"
    "../scipy/special/cdflib/cumfnc.f"
    "../scipy/special/cdflib/cumgam.f"
    "../scipy/special/cdflib/cumnbn.f"
    "../scipy/special/cdflib/cumnor.f"
    "../scipy/special/cdflib/cumpoi.f"
    "../scipy/special/cdflib/cumt.f"
    "../scipy/special/cdflib/cumtnc.f"
    "../scipy/special/cdflib/devlpl.f"
    "../scipy/special/cdflib/dinvnr.f"
    "../scipy/special/cdflib/dinvr.f"
    "../scipy/special/cdflib/dt1.f"
    "../scipy/special/cdflib/dzror.f"
    "../scipy/special/cdflib/erf.f"
    "../scipy/special/cdflib/erfc1.f"
    "../scipy/special/cdflib/esum.f"
    "../scipy/special/cdflib/exparg.f"
    "../scipy/special/cdflib/fpser.f"
    "../scipy/special/cdflib/gam1.f"
    "../scipy/special/cdflib/gaminv.f"
    "../scipy/special/cdflib/gamln.f"
    "../scipy/special/cdflib/gamln1.f"
    "../scipy/special/cdflib/gamma_fort.f"
    "../scipy/special/cdflib/grat1.f"
    "../scipy/special/cdflib/gratio.f"
    "../scipy/special/cdflib/gsumln.f"
    "../scipy/special/cdflib/ipmpar.f"
    "../scipy/special/cdflib/psi_fort.f"
    "../scipy/special/cdflib/rcomp.f"
    "../scipy/special/cdflib/rexp.f"
    "../scipy/special/cdflib/rlog.f"
    "../scipy/special/cdflib/rlog1.f"
    "../scipy/special/cdflib/spmpar.f"
    "../scipy/special/cdflib/stvaln.f"
)

SET(SOURCE_FILES_CDFLIB_RTLIB
    "../scipy/special/cdflib/algdiv.f"
)

SET(SOURCE_FILES_MINPACK2
    "../scipy/optimize/minpack2/dcsrch.f"
    "../scipy/optimize/minpack2/dcstep.f"
)

SET(SOURCE_FILES_MINPACK2_RTLIB
    "../scipy/optimize/minpack2/dcstep.f"
)

set(SOURCE_FILES_MACH_INTEGRATE
    "../scipy/integrate/mach/d1mach.f"
    "../scipy/integrate/mach/xerror.f"
)

SET(SOURCE_FILES_MACH_INTEGRATE_RTLIB
    "../scipy/integrate/mach/d1mach.f"
)

SET(SOURCE_FILES_QUADPACK
    "../scipy/integrate/quadpack/dqag.f"
    "../scipy/integrate/quadpack/dqage.f"
    "../scipy/integrate/quadpack/dqagi.f"
    "../scipy/integrate/quadpack/dqagie.f"
    "../scipy/integrate/quadpack/dqagp.f"
    "../scipy/integrate/quadpack/dqagpe.f"
    "../scipy/integrate/quadpack/dqags.f"
    "../scipy/integrate/quadpack/dqagse.f"
    "../scipy/integrate/quadpack/dqawc.f"
    "../scipy/integrate/quadpack/dqawce.f"
    "../scipy/integrate/quadpack/dqawf.f"
    "../scipy/integrate/quadpack/dqawfe.f"
    "../scipy/integrate/quadpack/dqawo.f"
    "../scipy/integrate/quadpack/dqawoe.f"
    "../scipy/integrate/quadpack/dqaws.f"
    "../scipy/integrate/quadpack/dqawse.f"
    "../scipy/integrate/quadpack/dqc25c.f"
    "../scipy/integrate/quadpack/dqc25f.f"
    "../scipy/integrate/quadpack/dqc25s.f"
    "../scipy/integrate/quadpack/dqcheb.f"
    "../scipy/integrate/quadpack/dqelg.f"
    "../scipy/integrate/quadpack/dqk15.f"
    "../scipy/integrate/quadpack/dqk15i.f"
    "../scipy/integrate/quadpack/dqk15w.f"
    "../scipy/integrate/quadpack/dqk21.f"
    "../scipy/integrate/quadpack/dqk31.f"
    "../scipy/integrate/quadpack/dqk41.f"
    "../scipy/integrate/quadpack/dqk51.f"
    "../scipy/integrate/quadpack/dqk61.f"
    "../scipy/integrate/quadpack/dqmomo.f"
    "../scipy/integrate/quadpack/dqng.f"
    "../scipy/integrate/quadpack/dqpsrt.f"
    "../scipy/integrate/quadpack/dqwgtc.f"
    "../scipy/integrate/quadpack/dqwgtf.f"
    "../scipy/integrate/quadpack/dqwgts.f"
)

SET(SOURCE_FILES_QUADPACK_RTLIB
    "../scipy/integrate/quadpack/dqwgts.f"
)

foreach(SRC ${SOURCE_FILES_MINPACK})
    set_source_files_properties(${SRC} PROPERTIES
        COMPILE_OPTIONS "--fixed-form;--implicit-typing;--implicit-interface;--generate-object-code;--use-loop-variable-after-loop;--all-mangling;--mangle-underscore;--bindc-mangling;--legacy-array-sections")
endforeach()

get_source_file_property(SOURCE_COMPILE_OPTIONS ${SOURCE_FILES_MINPACK_RTLIB} COMPILE_OPTIONS)
list(APPEND SOURCE_COMPILE_OPTIONS --rtlib)
set_source_files_properties(${SOURCE_FILES_MINPACK_RTLIB} PROPERTIES
    COMPILE_OPTIONS "${SOURCE_COMPILE_OPTIONS}")

foreach(SRC ${SOURCE_FILES_SPECFUN})
    set_source_files_properties(${SRC} PROPERTIES
        COMPILE_OPTIONS "--fixed-form;--implicit-typing;--implicit-interface;--generate-object-code;--use-loop-variable-after-loop;--all-mangling;--mangle-underscore;--bindc-mangling;--legacy-array-sections;--rtlib")
endforeach()

foreach(SRC ${SOURCE_FILES_FITPACK})
    set_source_files_properties(${SRC} PROPERTIES
        COMPILE_OPTIONS "--fixed-form;--implicit-typing;--implicit-interface;--generate-object-code;--use-loop-variable-after-loop;--all-mangling;--mangle-underscore;--bindc-mangling;--legacy-array-sections")
endforeach()

get_source_file_property(SOURCE_COMPILE_OPTIONS ${SOURCE_FILES_FITPACK_RTLIB} COMPILE_OPTIONS)
list(APPEND SOURCE_COMPILE_OPTIONS --rtlib)
set_source_files_properties(${SOURCE_FILES_FITPACK_RTLIB} PROPERTIES
    COMPILE_OPTIONS "${SOURCE_COMPILE_OPTIONS}")

foreach(SRC ${SOURCE_FILES_AMOS})
    set_source_files_properties(${SRC} PROPERTIES
        COMPILE_OPTIONS "--fixed-form;--implicit-typing;--implicit-interface;--generate-object-code;--use-loop-variable-after-loop;--all-mangling;--mangle-underscore;--bindc-mangling;--legacy-array-sections")
endforeach()

get_source_file_property(SOURCE_COMPILE_OPTIONS ${SOURCE_FILES_AMOS_RTLIB} COMPILE_OPTIONS)
list(APPEND SOURCE_COMPILE_OPTIONS --rtlib)
set_source_files_properties(${SOURCE_FILES_AMOS_RTLIB} PROPERTIES
    COMPILE_OPTIONS "${SOURCE_COMPILE_OPTIONS}")

foreach(SRC ${SOURCE_FILES_MACH})
    set_source_files_properties(${SRC} PROPERTIES
        COMPILE_OPTIONS "--fixed-form;--implicit-typing;--implicit-interface;--generate-object-code;--use-loop-variable-after-loop;--all-mangling;--mangle-underscore;--bindc-mangling;--legacy-array-sections")
endforeach()

foreach(SRC ${SOURCE_FILES_CDFLIB})
    set_source_files_properties(${SRC} PROPERTIES
        COMPILE_OPTIONS "--fixed-form;--implicit-typing;--implicit-interface;--generate-object-code;--use-loop-variable-after-loop;--all-mangling;--mangle-underscore;--bindc-mangling;--legacy-array-sections")
endforeach()

get_source_file_property(SOURCE_COMPILE_OPTIONS ${SOURCE_FILES_CDFLIB_RTLIB} COMPILE_OPTIONS)
list(APPEND SOURCE_COMPILE_OPTIONS --rtlib)
set_source_files_properties(${SOURCE_FILES_CDFLIB_RTLIB} PROPERTIES
    COMPILE_OPTIONS "${SOURCE_COMPILE_OPTIONS}")

foreach(SRC ${SOURCE_FILES_MINPACK2})
    set_source_files_properties(${SRC} PROPERTIES
        COMPILE_OPTIONS "--fixed-form;--implicit-typing;--implicit-interface;--generate-object-code;--use-loop-variable-after-loop;--all-mangling;--mangle-underscore;--bindc-mangling;--legacy-array-sections")
endforeach()

get_source_file_property(SOURCE_COMPILE_OPTIONS ${SOURCE_FILES_MINPACK2_RTLIB} COMPILE_OPTIONS)
list(APPEND SOURCE_COMPILE_OPTIONS --rtlib)
set_source_files_properties(${SOURCE_FILES_MINPACK2_RTLIB} PROPERTIES
    COMPILE_OPTIONS "${SOURCE_COMPILE_OPTIONS}")

    foreach(SRC ${SOURCE_FILES_MACH_INTEGRATE})
    set_source_files_properties(${SRC} PROPERTIES
        COMPILE_OPTIONS "--fixed-form;--implicit-typing;--implicit-interface;--generate-object-code;--use-loop-variable-after-loop;--all-mangling;--mangle-underscore;--bindc-mangling;--legacy-array-sections")
endforeach()

get_source_file_property(SOURCE_COMPILE_OPTIONS ${SOURCE_FILES_MACH_INTEGRATE_RTLIB} COMPILE_OPTIONS)
list(APPEND SOURCE_COMPILE_OPTIONS --rtlib)
set_source_files_properties(${SOURCE_FILES_MACH_INTEGRATE_RTLIB} PROPERTIES
    COMPILE_OPTIONS "${SOURCE_COMPILE_OPTIONS}")

foreach(SRC ${SOURCE_FILES_QUADPACK})
    set_source_files_properties(${SRC} PROPERTIES
        COMPILE_OPTIONS "--fixed-form;--implicit-typing;--implicit-interface;--generate-object-code;--use-loop-variable-after-loop;--all-mangling;--mangle-underscore;--bindc-mangling;--legacy-array-sections")
endforeach()

get_source_file_property(SOURCE_COMPILE_OPTIONS ${SOURCE_FILES_QUADPACK_RTLIB} COMPILE_OPTIONS)
list(APPEND SOURCE_COMPILE_OPTIONS --rtlib)
set_source_files_properties(${SOURCE_FILES_QUADPACK_RTLIB} PROPERTIES
    COMPILE_OPTIONS "${SOURCE_COMPILE_OPTIONS}")

# Create shared libraries for each component (Minpack, Fitpack, Specfun)
add_library(minpack SHARED ${SOURCE_FILES_MINPACK})
set_target_properties(minpack PROPERTIES POSITION_INDEPENDENT_CODE OFF)

add_library(specfun SHARED ${SOURCE_FILES_SPECFUN})
set_target_properties(specfun PROPERTIES POSITION_INDEPENDENT_CODE OFF)

add_library(fitpack SHARED ${SOURCE_FILES_FITPACK})
set_target_properties(fitpack PROPERTIES POSITION_INDEPENDENT_CODE OFF)

add_library(amos SHARED ${SOURCE_FILES_AMOS})
set_target_properties(amos PROPERTIES POSITION_INDEPENDENT_CODE OFF)

add_library(mach SHARED ${SOURCE_FILES_MACH})
set_target_properties(mach PROPERTIES POSITION_INDEPENDENT_CODE OFF)

add_library(cdflib SHARED ${SOURCE_FILES_CDFLIB})
set_target_properties(cdflib PROPERTIES POSITION_INDEPENDENT_CODE OFF)

add_library(minpack2 SHARED ${SOURCE_FILES_MINPACK2})
set_target_properties(minpack2 PROPERTIES POSITION_INDEPENDENT_CODE OFF)

add_library(mach_integrate SHARED ${SOURCE_FILES_MACH_INTEGRATE})
set_target_properties(mach_integrate PROPERTIES POSITION_INDEPENDENT_CODE OFF)

add_library(quadpack SHARED ${SOURCE_FILES_QUADPACK})
set_target_properties(quadpack PROPERTIES POSITION_INDEPENDENT_CODE OFF)

# Link the libraries
message("${LFORTRAN_RUNTIME_LIBRARY_PATH}")
find_library(LFORTRAN_RUNTIME_LIBRARY lfortran_runtime PATH "${LFORTRAN_RUNTIME_LIBRARY_PATH}")
message("${LFORTRAN_RUNTIME_LIBRARY}")
add_library(lfortran_runtime SHARED IMPORTED)
set_property(TARGET lfortran_runtime PROPERTY IMPORTED_LOCATION "${LFORTRAN_RUNTIME_LIBRARY}")
target_link_libraries(minpack PRIVATE lfortran_runtime)
target_link_libraries(specfun PRIVATE lfortran_runtime)
target_link_libraries(fitpack PRIVATE lfortran_runtime)
target_link_libraries(amos PRIVATE lfortran_runtime)
target_link_libraries(mach PRIVATE lfortran_runtime)
target_link_libraries(cdflib PRIVATE lfortran_runtime)
target_link_libraries(minpack2 PRIVATE lfortran_runtime)
target_link_libraries(mach_integrate PRIVATE lfortran_runtime)
target_link_libraries(quadpack PRIVATE lfortran_runtime)
# Set C as the link language
set_target_properties(minpack specfun fitpack amos mach cdflib minpack2 mach_integrate quadpack PROPERTIES LINKER_LANGUAGE C)

# Rename the shared library to have a .so extension
set_target_properties(minpack PROPERTIES OUTPUT_NAME minpack SUFFIX ".so")
set_target_properties(specfun PROPERTIES OUTPUT_NAME specfun SUFFIX ".so")
set_target_properties(fitpack PROPERTIES OUTPUT_NAME fitpack SUFFIX ".so")
set_target_properties(amos PROPERTIES OUTPUT_NAME amos SUFFIX ".so")
set_target_properties(mach PROPERTIES OUTPUT_NAME mach SUFFIX ".so")
set_target_properties(cdflib PROPERTIES OUTPUT_NAME cdflib SUFFIX ".so")
set_target_properties(minpack2 PROPERTIES OUTPUT_NAME minpack2 SUFFIX ".so")
set_target_properties(mach_integrate PROPERTIES OUTPUT_NAME mach_integrate SUFFIX ".so")
set_target_properties(quadpack PROPERTIES OUTPUT_NAME quadpack SUFFIX ".so")

# Specify installation paths
install(TARGETS minpack DESTINATION ${MINPACK_DIRECTORY})
install(TARGETS minpack DESTINATION $ENV{CONDA_PREFIX}/lib)

install(TARGETS specfun DESTINATION ${SPECFUN_DIRECTORY})
install(TARGETS specfun DESTINATION $ENV{CONDA_PREFIX}/lib)

install(TARGETS fitpack DESTINATION ${FITPACK_DIRECTORY})
install(TARGETS fitpack DESTINATION $ENV{CONDA_PREFIX}/lib)

install(TARGETS amos DESTINATION ${SPECFUN_DIRECTORY})
install(TARGETS amos DESTINATION $ENV{CONDA_PREFIX}/lib)

install(TARGETS mach DESTINATION ${SPECFUN_DIRECTORY})
install(TARGETS mach DESTINATION $ENV{CONDA_PREFIX}/lib)

install(TARGETS cdflib DESTINATION ${SPECFUN_DIRECTORY})
install(TARGETS cdflib DESTINATION $ENV{CONDA_PREFIX}/lib)

install(TARGETS minpack2 DESTINATION ${MINPACK_DIRECTORY})
install(TARGETS minpack2 DESTINATION $ENV{CONDA_PREFIX}/lib)

install(TARGETS mach_integrate DESTINATION ${INTEGRATE_DIRECTORY})
install(TARGETS mach_integrate DESTINATION $ENV{CONDA_PREFIX}/lib)

install(TARGETS quadpack DESTINATION ${INTEGRATE_DIRECTORY})
install(TARGETS quadpack DESTINATION $ENV{CONDA_PREFIX}/lib)

