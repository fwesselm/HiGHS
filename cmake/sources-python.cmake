set(include_dirs_python
    ${CMAKE_SOURCE_DIR}/extern
    ${CMAKE_SOURCE_DIR}/extern/filereader
    ${CMAKE_SOURCE_DIR}/extern/pdqsort
    ${CMAKE_SOURCE_DIR}/extern/zstr
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/interfaces
    ${CMAKE_SOURCE_DIR}/src/io
    ${CMAKE_SOURCE_DIR}/src/ipm
    ${CMAKE_SOURCE_DIR}/src/ipm/ipx
    ${CMAKE_SOURCE_DIR}/src/ipm/basiclu
    ${CMAKE_SOURCE_DIR}/src/lp_data
    ${CMAKE_SOURCE_DIR}/src/mip
    ${CMAKE_SOURCE_DIR}/src/model
    ${CMAKE_SOURCE_DIR}/src/parallel
    ${CMAKE_SOURCE_DIR}/src/pdlp
    ${CMAKE_SOURCE_DIR}/src/pdlp/cupdlp
    ${CMAKE_SOURCE_DIR}/src/presolve
    ${CMAKE_SOURCE_DIR}/src/qpsolver
    ${CMAKE_SOURCE_DIR}/src/simplex
    ${CMAKE_SOURCE_DIR}/src/test
    ${CMAKE_SOURCE_DIR}/src/util
    $<BUILD_INTERFACE:${HIGHS_BINARY_DIR}>)

set(cupdlp_sources_python
  src/pdlp/cupdlp/cupdlp_cs.c
  src/pdlp/cupdlp/cupdlp_linalg.c
  src/pdlp/cupdlp/cupdlp_proj.c
  src/pdlp/cupdlp/cupdlp_restart.c
  src/pdlp/cupdlp/cupdlp_scaling.c
  src/pdlp/cupdlp/cupdlp_solver.c
  src/pdlp/cupdlp/cupdlp_step.c
  src/pdlp/cupdlp/cupdlp_utils.c)

set(cupdlp_headers_python
  src/pdlp/cupdlp/cupdlp_cs.h
  src/pdlp/cupdlp/cupdlp_defs.h
  src/pdlp/cupdlp/cupdlp_linalg.h
  src/pdlp/cupdlp/cupdlp_proj.h
  src/pdlp/cupdlp/cupdlp_restart.h
  src/pdlp/cupdlp/cupdlp_scaling.h
  src/pdlp/cupdlp/cupdlp_solver.h
  src/pdlp/cupdlp/cupdlp_step.h
  src/pdlp/cupdlp/cupdlp_utils.c)

set(cuda_sources_python
  pdlp/cupdlp/cuda/cupdlp_cuda_kernels.cu
  pdlp/cupdlp/cuda/cupdlp_cuda_kernels.cuh
  pdlp/cupdlp/cuda/cupdlp_cudalinalg.cuh
  pdlp/cupdlp/cuda/cupdlp_cudalinalg.cu)

set(basiclu_sources_python
  src/ipm/basiclu/basiclu_factorize.c
  src/ipm/basiclu/basiclu_get_factors.c
  src/ipm/basiclu/basiclu_initialize.c
  src/ipm/basiclu/basiclu_object.c
  src/ipm/basiclu/basiclu_solve_dense.c
  src/ipm/basiclu/basiclu_solve_for_update.c
  src/ipm/basiclu/basiclu_solve_sparse.c
  src/ipm/basiclu/basiclu_update.c
  src/ipm/basiclu/lu_build_factors.c
  src/ipm/basiclu/lu_condest.c
  src/ipm/basiclu/lu_dfs.c
  src/ipm/basiclu/lu_factorize_bump.c
  src/ipm/basiclu/lu_file.c
  src/ipm/basiclu/lu_garbage_perm.c
  src/ipm/basiclu/lu_initialize.c
  src/ipm/basiclu/lu_internal.c
  src/ipm/basiclu/lu_markowitz.c
  src/ipm/basiclu/lu_matrix_norm.c
  src/ipm/basiclu/lu_pivot.c
  src/ipm/basiclu/lu_residual_test.c
  src/ipm/basiclu/lu_setup_bump.c
  src/ipm/basiclu/lu_singletons.c
  src/ipm/basiclu/lu_solve_dense.c
  src/ipm/basiclu/lu_solve_for_update.c
  src/ipm/basiclu/lu_solve_sparse.c
  src/ipm/basiclu/lu_solve_symbolic.c
  src/ipm/basiclu/lu_solve_triangular.c
  src/ipm/basiclu/lu_update.c)

set(basiclu_headers_python
  src/ipm/basiclu/basiclu_factorize.h
  src/ipm/basiclu/basiclu_get_factors.h
  src/ipm/basiclu/basiclu_initialize.h
  src/ipm/basiclu/basiclu_obj_factorize.h
  src/ipm/basiclu/basiclu_obj_free.h
  src/ipm/basiclu/basiclu_obj_get_factors.h
  src/ipm/basiclu/basiclu_obj_initialize.h
  src/ipm/basiclu/basiclu_obj_solve_dense.h
  src/ipm/basiclu/basiclu_obj_solve_for_update.h
  src/ipm/basiclu/basiclu_obj_solve_sparse.h
  src/ipm/basiclu/basiclu_obj_update.h
  src/ipm/basiclu/basiclu_object.h
  src/ipm/basiclu/basiclu_solve_dense.h
  src/ipm/basiclu/basiclu_solve_for_update.h
  src/ipm/basiclu/basiclu_solve_sparse.h
  src/ipm/basiclu/basiclu_update.h
  src/ipm/basiclu/basiclu.h
  src/ipm/basiclu/lu_def.h
  src/ipm/basiclu/lu_file.h
  src/ipm/basiclu/lu_internal.h
  src/ipm/basiclu/lu_list.h)

set(ipx_sources_python
  src/ipm/ipx/basiclu_kernel.cc
  src/ipm/ipx/basiclu_wrapper.cc
  src/ipm/ipx/basis.cc
  src/ipm/ipx/conjugate_residuals.cc
  src/ipm/ipx/control.cc
  src/ipm/ipx/crossover.cc
  src/ipm/ipx/diagonal_precond.cc
  src/ipm/ipx/forrest_tomlin.cc
  src/ipm/ipx/guess_basis.cc
  src/ipm/ipx/indexed_vector.cc
  src/ipm/ipx/info.cc
  src/ipm/ipx/ipm.cc
  src/ipm/ipx/ipx_c.cc
  src/ipm/ipx/iterate.cc
  src/ipm/ipx/kkt_solver_basis.cc
  src/ipm/ipx/kkt_solver_diag.cc
  src/ipm/ipx/kkt_solver.cc
  src/ipm/ipx/linear_operator.cc
  src/ipm/ipx/lp_solver.cc
  src/ipm/ipx/lu_factorization.cc
  src/ipm/ipx/lu_update.cc
  src/ipm/ipx/maxvolume.cc
  src/ipm/ipx/model.cc
  src/ipm/ipx/normal_matrix.cc
  src/ipm/ipx/sparse_matrix.cc
  src/ipm/ipx/sparse_utils.cc
  src/ipm/ipx/splitted_normal_matrix.cc
  src/ipm/ipx/starting_basis.cc
  src/ipm/ipx/symbolic_invert.cc
  src/ipm/ipx/timer.cc
  src/ipm/ipx/utils.cc)
  
  set(ipx_headers_python
  src/ipm/ipx/basiclu_kernel.h
  src/ipm/ipx/basiclu_wrapper.h
  src/ipm/ipx/basis.h
  src/ipm/ipx/conjugate_residuals.h
  src/ipm/ipx/control.h
  src/ipm/ipx/crossover.h
  src/ipm/ipx/diagonal_precond.h
  src/ipm/ipx/forrest_tomlin.h
  src/ipm/ipx/guess_basis.h
  src/ipm/ipx/indexed_vector.h
  src/ipm/ipx/info.h
  src/ipm/ipx/ipm.h
  src/ipm/ipx/ipx_c.h
  src/ipm/ipx/ipx_config.h
  src/ipm/ipx/ipx_info.h
  src/ipm/ipx/ipx_internal.h
  src/ipm/ipx/ipx_parameters.h
  src/ipm/ipx/ipx_status.h
  src/ipm/ipx/iterate.h
  src/ipm/ipx/kkt_solver_basis.h
  src/ipm/ipx/kkt_solver_diag.h
  src/ipm/ipx/kkt_solver.h
  src/ipm/ipx/linear_operator.h
  src/ipm/ipx/lp_solver.h
  src/ipm/ipx/lu_factorization.h
  src/ipm/ipx/lu_update.h
  src/ipm/ipx/maxvolume.h
  src/ipm/ipx/model.h
  src/ipm/ipx/multistream.h
  src/ipm/ipx/normal_matrix.h
  src/ipm/ipx/power_method.h
  src/ipm/ipx/sparse_matrix.h
  src/ipm/ipx/sparse_utils.h
  src/ipm/ipx/splitted_normal_matrix.h
  src/ipm/ipx/starting_basis.h
  src/ipm/ipx/symbolic_invert.h
  src/ipm/ipx/timer.h
  src/ipm/ipx/utils.h)

set(highs_sources_python
    extern/filereaderlp/reader.cpp
    src/interfaces/highs_c_api.cpp
    src/io/Filereader.cpp
    src/io/FilereaderEms.cpp
    src/io/FilereaderLp.cpp
    src/io/FilereaderMps.cpp
    src/io/HighsIO.cpp
    src/io/HMpsFF.cpp
    src/io/HMPSIO.cpp
    src/io/LoadOptions.cpp
    src/ipm/IpxWrapper.cpp
    src/lp_data/Highs.cpp
    src/lp_data/HighsCallback.cpp
    src/lp_data/HighsDebug.cpp
    src/lp_data/HighsIis.cpp
    src/lp_data/HighsInfo.cpp
    src/lp_data/HighsInfoDebug.cpp
    src/lp_data/HighsInterface.cpp
    src/lp_data/HighsLp.cpp
    src/lp_data/HighsLpUtils.cpp
    src/lp_data/HighsModelUtils.cpp
    src/lp_data/HighsOptions.cpp
    src/lp_data/HighsRanging.cpp
    src/lp_data/HighsSolution.cpp
    src/lp_data/HighsSolutionDebug.cpp
    src/lp_data/HighsSolve.cpp
    src/lp_data/HighsStatus.cpp
    src/mip/HighsCliqueTable.cpp
    src/mip/HighsConflictPool.cpp
    src/mip/HighsCutGeneration.cpp
    src/mip/HighsCutPool.cpp
    src/mip/HighsDebugSol.cpp
    src/mip/HighsDomain.cpp
    src/mip/HighsDynamicRowMatrix.cpp
    src/mip/HighsGFkSolve.cpp
    src/mip/HighsImplications.cpp
    src/mip/HighsLpAggregator.cpp
    src/mip/HighsLpRelaxation.cpp
    src/mip/HighsMipAnalysis.cpp
    src/mip/HighsMipSolver.cpp
    src/mip/HighsMipSolverData.cpp
    src/mip/HighsModkSeparator.cpp
    src/mip/HighsNodeQueue.cpp
    src/mip/HighsObjectiveFunction.cpp
    src/mip/HighsPathSeparator.cpp
    src/mip/HighsPrimalHeuristics.cpp
    src/mip/HighsPseudocost.cpp
    src/mip/HighsRedcostFixing.cpp
    src/mip/HighsSearch.cpp
    src/mip/HighsSeparation.cpp
    src/mip/HighsSeparator.cpp
    src/mip/HighsTableauSeparator.cpp
    src/mip/HighsTransformedLp.cpp
    src/model/HighsHessian.cpp
    src/model/HighsHessianUtils.cpp
    src/model/HighsModel.cpp
    src/parallel/HighsTaskExecutor.cpp
    src/pdlp/CupdlpWrapper.cpp
    src/presolve/HighsPostsolveStack.cpp
    src/presolve/HighsSymmetry.cpp
    src/presolve/HPresolve.cpp
    src/presolve/HPresolveAnalysis.cpp
    src/presolve/ICrash.cpp
    src/presolve/ICrashUtil.cpp
    src/presolve/ICrashX.cpp
    src/presolve/PresolveComponent.cpp
    src/qpsolver/a_asm.cpp
    src/qpsolver/a_quass.cpp
    src/qpsolver/basis.cpp
    src/qpsolver/perturbation.cpp
    src/qpsolver/quass.cpp
    src/qpsolver/ratiotest.cpp
    src/qpsolver/scaling.cpp
    src/simplex/HEkk.cpp
    src/simplex/HEkkControl.cpp
    src/simplex/HEkkDebug.cpp
    src/simplex/HEkkDual.cpp
    src/simplex/HEkkDualMulti.cpp
    src/simplex/HEkkDualRHS.cpp
    src/simplex/HEkkDualRow.cpp
    src/simplex/HEkkInterface.cpp
    src/simplex/HEkkPrimal.cpp
    src/simplex/HighsSimplexAnalysis.cpp
    src/simplex/HSimplex.cpp
    src/simplex/HSimplexDebug.cpp
    src/simplex/HSimplexNla.cpp
    src/simplex/HSimplexNlaDebug.cpp
    src/simplex/HSimplexNlaFreeze.cpp
    src/simplex/HSimplexNlaProductForm.cpp
    src/simplex/HSimplexReport.cpp
    src/test/KktCh2.cpp
    src/test/DevKkt.cpp
    src/util/HFactor.cpp
    src/util/HFactorDebug.cpp
    src/util/HFactorExtend.cpp
    src/util/HFactorRefactor.cpp
    src/util/HFactorUtils.cpp
    src/util/HighsHash.cpp
    src/util/HighsLinearSumBounds.cpp
    src/util/HighsMatrixPic.cpp
    src/util/HighsMatrixUtils.cpp
    src/util/HighsSort.cpp
    src/util/HighsSparseMatrix.cpp
    src/util/HighsUtils.cpp
    src/util/HSet.cpp
    src/util/HVectorBase.cpp
    src/util/stringutil.cpp)

set(highs_headers_python
    extern/filereaderlp/builder.hpp
    extern/filereaderlp/def.hpp
    extern/filereaderlp/model.hpp
    extern/filereaderlp/reader.hpp
    extern/pdqsort/pdqsort.h
    src/interfaces/highs_c_api.h
    src/io/Filereader.h
    src/io/FilereaderEms.h
    src/io/FilereaderLp.h
    src/io/FilereaderMps.h
    src/io/HighsIO.h
    src/io/HMpsFF.h
    src/io/HMPSIO.h
    src/io/LoadOptions.h
    src/ipm/IpxSolution.h
    src/ipm/IpxWrapper.h
    src/lp_data/HConst.h
    src/lp_data/HighsAnalysis.h
    src/lp_data/HighsCallback.h
    src/lp_data/HighsCallbackStruct.h
    src/lp_data/HighsDebug.h
    src/lp_data/HighsIis.h
    src/lp_data/HighsInfo.h
    src/lp_data/HighsInfoDebug.h
    src/lp_data/HighsLp.h
    src/lp_data/HighsLpSolverObject.h
    src/lp_data/HighsLpUtils.h
    src/lp_data/HighsModelUtils.h
    src/lp_data/HighsOptions.h
    src/lp_data/HighsRanging.h
    src/lp_data/HighsSolution.h
    src/lp_data/HighsSolutionDebug.h
    src/lp_data/HighsSolve.h
    src/lp_data/HighsStatus.h
    src/lp_data/HStruct.h
    src/mip/HighsCliqueTable.h
    src/mip/HighsConflictPool.h
    src/mip/HighsCutGeneration.h
    src/mip/HighsCutPool.h
    src/mip/HighsDebugSol.h
    src/mip/HighsDomain.h
    src/mip/HighsDomainChange.h
    src/mip/HighsDynamicRowMatrix.h
    src/mip/HighsGFkSolve.h
    src/mip/HighsImplications.h
    src/mip/HighsLpAggregator.h
    src/mip/HighsLpRelaxation.h
    src/mip/HighsMipAnalysis.h
    src/mip/HighsMipSolver.h
    src/mip/HighsMipSolverData.h
    src/mip/HighsModkSeparator.h
    src/mip/HighsNodeQueue.h
    src/mip/HighsObjectiveFunction.h
    src/mip/HighsPathSeparator.h
    src/mip/HighsPrimalHeuristics.h
    src/mip/HighsPseudocost.h
    src/mip/HighsRedcostFixing.h
    src/mip/HighsSearch.h
    src/mip/HighsSeparation.h
    src/mip/HighsSeparator.h
    src/mip/HighsTableauSeparator.h
    src/mip/HighsTransformedLp.h
    src/mip/MipTimer.h
    src/model/HighsHessian.h
    src/model/HighsHessianUtils.h
    src/model/HighsModel.h
    src/parallel/HighsBinarySemaphore.h
    src/parallel/HighsCacheAlign.h
    src/parallel/HighsCombinable.h
    src/parallel/HighsMutex.h
    src/parallel/HighsParallel.h
    src/parallel/HighsRaceTimer.h
    src/parallel/HighsSchedulerConstants.h
    src/parallel/HighsSpinMutex.h
    src/parallel/HighsSplitDeque.h
    src/parallel/HighsTask.h
    src/parallel/HighsTaskExecutor.h
    src/pdlp/CupdlpWrapper.h
    src/presolve/HighsPostsolveStack.h
    src/presolve/HighsSymmetry.h
    src/presolve/HPresolve.h
    src/presolve/HPresolveAnalysis.h
    src/presolve/ICrash.h
    src/presolve/ICrashUtil.h
    src/presolve/ICrashX.h
    src/presolve/PresolveComponent.h
    src/qpsolver/a_asm.hpp
    src/qpsolver/a_quass.hpp
    src/qpsolver/basis.hpp
    src/qpsolver/crashsolution.hpp
    src/qpsolver/dantzigpricing.hpp
    src/qpsolver/devexpricing.hpp
    src/qpsolver/eventhandler.hpp
    src/qpsolver/factor.hpp
    src/qpsolver/feasibility_bounded.hpp
    src/qpsolver/feasibility_highs.hpp
    src/qpsolver/gradient.hpp
    src/qpsolver/instance.hpp
    src/qpsolver/matrix.hpp
    src/qpsolver/perturbation.hpp
    src/qpsolver/pricing.hpp
    src/qpsolver/qpconst.hpp
    src/qpsolver/qpvector.hpp
    src/qpsolver/quass.hpp
    src/qpsolver/ratiotest.hpp
    src/qpsolver/runtime.hpp
    src/qpsolver/scaling.hpp
    src/qpsolver/settings.hpp
    src/qpsolver/snippets.hpp
    src/qpsolver/statistics.hpp
    src/qpsolver/steepestedgepricing.hpp
    src/simplex/HApp.h
    src/simplex/HEkk.h
    src/simplex/HEkkDual.h
    src/simplex/HEkkDualRHS.h
    src/simplex/HEkkDualRow.h
    src/simplex/HEkkPrimal.h
    src/simplex/HighsSimplexAnalysis.h
    src/simplex/HSimplex.h
    src/simplex/HSimplexDebug.h
    src/simplex/HSimplexNla.h
    src/simplex/HSimplexReport.h
    src/simplex/SimplexConst.h
    src/simplex/SimplexStruct.h
    src/simplex/SimplexTimer.h
    src/test/DevKkt.h
    src/test/KktCh2.h
    src/util/FactorTimer.h
    src/util/HFactor.h
    src/util/HFactorConst.h
    src/util/HFactorDebug.h
    src/util/HighsCDouble.h
    src/util/HighsComponent.h
    src/util/HighsDataStack.h
    src/util/HighsDisjointSets.h
    src/util/HighsHash.h
    src/util/HighsHashTree.h
    src/util/HighsInt.h
    src/util/HighsIntegers.h
    src/util/HighsLinearSumBounds.h
    src/util/HighsMatrixPic.h
    src/util/HighsMatrixSlice.h
    src/util/HighsMatrixUtils.h
    src/util/HighsMemoryAllocation.h
    src/util/HighsRandom.h
    src/util/HighsRbTree.h
    src/util/HighsSort.h
    src/util/HighsSparseMatrix.h
    src/util/HighsSparseVectorSum.h
    src/util/HighsSplay.h
    src/util/HighsTimer.h
    src/util/HighsUtils.h
    src/util/HSet.h
    src/util/HVector.h
    src/util/HVectorBase.h
    src/util/stringutil.h
    src/Highs.h
  )