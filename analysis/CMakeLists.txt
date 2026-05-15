project(analysis)

# Analysis code was intentionally archived; keep this subdirectory
# valid for top-level CMake configuration.

set(
  analysis_SCRIPTS
  plot_trigger_timing.py
  compare_timing_resolution.py
  plot_cfd_per_trigger_interp.py
)

foreach(_script ${analysis_SCRIPTS})
  if(EXISTS "${PROJECT_SOURCE_DIR}/${_script}")
    configure_file(
      ${PROJECT_SOURCE_DIR}/${_script}
      ${PROJECT_BINARY_DIR}/${_script}
      COPYONLY
    )
  endif()
endforeach()
