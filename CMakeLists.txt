# Build using:
# FC=lfortran cmake -DCMAKE_Fortran_FLAGS=--verbose .
cmake_minimum_required(VERSION 3.10)

project(LFortranScipyBuild)
enable_language(Fortran)

set(SPECFUN_DIRECTORY ${CMAKE_BINARY_DIR}/../scipy/special/)
set(MINPACK_DIRECTORY ${CMAKE_BINARY_DIR}/../scipy/optimize/)
set(FITPACK_DIRECTORY ${CMAKE_BINARY_DIR}/../scipy/interpolate/)



# Specify the source files explicitly for each component Specfun, Minpack, Fitpack
set(SOURCE_FILES_SPECFUN
    "../scipy/special/specfun/specfun.f"
    # Add more source files here if needed
)

set(SOURCE_FILES_MINPACK
    "scipy/optimize/minpack/dogleg.f"
    "scipy/optimize/minpack/dpmpar.f"
    "scipy/optimize/minpack/enorm.f"
    "scipy/optimize/minpack/fdjac1.f"
    "scipy/optimize/minpack/fdjac2.f"
    "scipy/optimize/minpack/hybrd.f"
    "scipy/optimize/minpack/hybrd1.f"
    "scipy/optimize/minpack/hybrj.f"
    "scipy/optimize/minpack/hybrj1.f"
    "scipy/optimize/minpack/lmder.f"
    "scipy/optimize/minpack/lmder1.f"
    "scipy/optimize/minpack/lmdif.f"
    "scipy/optimize/minpack/lmdif1.f"
    "scipy/optimize/minpack/lmpar.f"
    "scipy/optimize/minpack/lmstr.f"
    "scipy/optimize/minpack/lmstr1.f"
    "scipy/optimize/minpack/qform.f"
    "scipy/optimize/minpack/qrfac.f"
    "scipy/optimize/minpack/qrsolv.f"
    "scipy/optimize/minpack/r1mpyq.f"
    "scipy/optimize/minpack/r1updt.f"
    "scipy/optimize/minpack/rwupdt.f"
    "scipy/optimize/minpack/chkder.f"
)

set(SOURCE_FILES_MINPACK_RTLIB
    "scipy/optimize/minpack/chkder.f"
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
)

SET(SOURCE_FILES_FITPACK_RTLIB
    "../scipy/interpolate/fitpack/evapol.f"
)

foreach(SRC ${SOURCE_FILES_MINPACK})
    set_source_files_properties(${SRC} PROPERTIES
        COMPILE_OPTIONS "--fixed-form;--implicit-typing;--implicit-interface;--generate-object-code;--use-loop-variable-after-loop;--all-mangling;--mangle-underscore;--bindc-mangling")
endforeach()

get_source_file_property(SOURCE_COMPILE_OPTIONS ${SOURCE_FILES_MINPACK_RTLIB} COMPILE_OPTIONS)
list(APPEND SOURCE_COMPILE_OPTIONS --rtlib)
message("XX: ${SOURCE_COMPILE_OPTIONS}")
set_source_files_properties(${SOURCE_FILES_MINPACK_RTLIB} PROPERTIES
    COMPILE_OPTIONS "${SOURCE_COMPILE_OPTIONS}")

# Create shared libraries for each component (Minpack, Fitpack, Specfun)
add_library(minpack SHARED ${SOURCE_FILES_MINPACK})
set_target_properties(minpack PROPERTIES POSITION_INDEPENDENT_CODE OFF)


# Link the libraries
target_link_libraries(minpack PRIVATE ${LFORTRAN_RUNTIME_LIBRARY_PATH})

# Set Fortran as the link language
set_target_properties(minpack PROPERTIES LINKER_LANGUAGE Fortran)

# Rename the shared library to have a .so extension
set_target_properties(minpack PROPERTIES OUTPUT_NAME minpack SUFFIX ".so")

# Specify installation paths
install(TARGETS minpack DESTINATION ${MINPACK_DIRECTORY})
install(TARGETS minpack DESTINATION ${CONDA_PREFIX_LIB_DIR})
