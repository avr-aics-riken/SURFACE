-- Build unit test
newoption {
   trigger     = "with-unittest",
   description = "Build unit test."
}

-- Enable profle setting for K/FX10
newoption {
   trigger     = "with-k-profile",
   description = "Insert profile flag for K/FX10."
}


-- Enable autoamtic screen parallel(valid with --with-mpi flag)
newoption {
   trigger     = "with-screen-parallel",
   description = "Enable automatic screen parallel rendering(also must specify --with-mpi)."
}

newoption {
   trigger     = "with-debugtrace",
   description = "Enable debug trace for better API debugging in runtime."
}

-- OpenMP
newoption {
   trigger     = "with-openmp",
   description = "Use OpenMP."
}

-- MPI
newoption {
   trigger     = "with-mpi",
   description = "Use MPI."
}

-- Tentative
newoption {
   trigger     = "with-glsl-opt",
   description = "Use Optimized GLSL codepath."
}

--
-- Defines custom platform.
--

function newplatform(plf)
    local name = plf.name
    local description = plf.description
 
    -- Register new platform
    premake.platforms[name] = {
        cfgsuffix = "_"..name,
        iscrosscompiler = true
    }
 
    if premake.option.list["platform"] == nil then
      -- Windows platform? Do nothing.
      return
    end

    -- Allow use of new platform in --platfroms
    table.insert(premake.option.list["platform"].allowed, { name, description })
    table.insert(premake.fields.platforms.allowed, name)
 
    -- Add compiler support
    -- gcc
    premake.gcc.platforms[name] = plf.gcc
    --other compilers (?)
end

newplatform {
    name = "k-cross",
    description = "K/FX10 cross compiler",
    gcc = {
        cc = "fccpx",
        cxx = "FCCpx",
        cppflags = "-MMD -Xg"
    }
}

newplatform {
    name = "k-cross-mpi",
    description = "K/FX10 MPI cross compiler",
    gcc = {
        cc = "mpifccpx",
        cxx = "mpiFCCpx",
        cppflags = "-MMD -Xg"
    }
}

newplatform {
    name = "k-native",
    description = "K/FX10 native compiler",
    gcc = {
        cc = "fcc",
        cxx = "FCC",
        cppflags = "-MMD -Xg"
    }
}

--
-- End platform definition.
--

render_sources = {
   "../render/accel_triangle.cc",
   "../render/accel_line.cc",
   "../render/accel_particle.cc",
   "../render/accel_tetra.cc",
   "../render/accel_volume.cc",
   "../render/camera.cc",
   "../render/prim_mesh.cc",
   "../render/texture.cc",
   "../render/tinymt64.cpp",
   "../render/toplevel_bvh.cc",
   -- "../render/logutil.cc",
   "../render/prefix_tree_util.cc",
   }

gles_sources = {
   --"gles_c_api.cc",
   "gles_common.cc",
   "gles_common.h",
   "gles_context.cc",
   "gles_context.h",
   "gles_buffer.cc",
   "gles_fragment_op.cc",
   "gles_graphics_state.cc",
   "gles_handle_allocator.h",
   --"gles_log.cc",
   "gles_accel_builder.h",
   "gles_accel_builder.cc",
   "gles_primitive.cc",
   "gles_program.cc",
   "gles_raytrace_engine.cc",
   "gles_resource_manager.h",
   "gles_resource_manager.cc",
   "gles_shader.cc",
   "gles_special.cc",
   "gles_texture.cc",
   "gles_vertex.cc",
   "gles_glsl_runtime.cc",
   "gles_render_graph.cc",
   }

gtest_sources = {
   "../third_party/gtest-1.7.0/src/gtest-all.cc",
}

test_sources = {
   "../test/cctest/test-dummy.cc",
   "../test/cctest/test-buffer.cc",
   "../test/cctest/test-graphics_state.cc",
   "../test/cctest/test-lsgl_ext.cc",
   "../test/cctest/test-shader.cc",
   "../test/cctest/test-draw.cc",
   "../test/cctest/test-vertexattrib.cc",
   "../test/cctest/test-matrix.cc",
   "../test/cctest/test-simdclass.cc",
   "../test/cctest/test-texture.cc",
   "../test/cctest/test-accel-particle.cc",
   "../test/cctest/test-prefix-tree.cc",
   "../test/cctest/test-main.cc",
}

-- premake4.lua
solution "libLSGLESSolution"
   configurations { "Release", "Debug" }

   k_platform = false

   if (os.is("windows")) then

      platforms { "x32", "x64", "native" }
      files { "gles_c_api.cc" }

   else

      -- default to x64 platform if none has been specified
      _OPTIONS.platform = _OPTIONS.platform or 'x64'

      -- overwrite the native platform with the options::platform
      premake.gcc.platforms['Native'] = premake.gcc.platforms[_OPTIONS.platform]

      k_platform = (_OPTIONS.platform == 'k-cross' or _OPTIONS.platform == 'k-native' or _OPTIONS.platform == 'k-cross-mpi')

      platforms { "native", "x32", "x64" }
   end

   -- A project defines one build target
   project "LSGLES"
      kind "SharedLib"
      language "C++"

      files { gles_sources, render_sources }

      includedirs {
         "../include"
      }

      if _OPTIONS['with-debugtrace'] then
         defines { "LSGL_DEBUG_TRACE" }
      end

      -- MacOSX. Guess we use gcc.
      configuration { "macosx", "gmake" }

         if _OPTIONS['with-glsl-opt'] then
            defines { 'LSGL_OPTIMIZE_GLSL' }
         end

         -- Screen parallel
         if _OPTIONS['with-screen-parallel'] then
            defines { "LSGL_ENABLE_SCREEN_PARALLEL" }
         end

         -- MPI
         if _OPTIONS['with-mpi'] then
            defines { "LSGL_ENABLE_MPI" }
         end

         -- gcc openmp
         if _OPTIONS['with-openmp'] then
            buildoptions { "-fopenmp" }
            --linkoptions { "-fopenmp" }

            -- statically link libgomp.
            linkoptions { "`${CXX} -print-file-name=libgomp.a`" }

            -- statically link libgcc and libstdc++
            --linkoptions { "-static-libgcc", "-static-libstdc++" }
         end


      -- Windows specific. Assumes visual studio.
      configuration { "windows" }

         -- Add MemoryModule code.
         files { '../deps/MemoryModule/MemoryModule.c' }
         includedirs { '../deps/MemoryModule' }

         defines { "LSGL_EXPORT" }

         includedirs { './compat' } -- stdint.h 
         defines { 'NOMINMAX', '_LARGEFILE_SOURCE', '_FILE_OFFSET_BITS=64' }
         defines { '__STDC_CONSTANT_MACROS', '__STDC_LIMIT_MACROS' } -- c99

         -- Screen parallel
         if _OPTIONS['with-screen-parallel'] then
            defines { "LSGL_ENABLE_SCREEN_PARALLEL" }
         end

         -- MPI
         if _OPTIONS['with-mpi'] then
            defines { "LSGL_ENABLE_MPI" }
         end

      -- Linux specific
      configuration { "linux" }

         if _OPTIONS['with-glsl-opt'] then
            defines { 'LSGL_OPTIMIZE_GLSL' }
            if k_platform then 
               buildoptions { '-x100', '-Kilfunc', '-Klib' }
            end
         end

         if _OPTIONS['with-k-profile'] then
            defines { "LSGL_ENABLE_K_PROFILE" }
         end

         defines { '__STDC_CONSTANT_MACROS', '__STDC_LIMIT_MACROS' } -- c99
         defines { '_LARGEFILE_SOURCE', '_FILE_OFFSET_BITS=64' }

         -- Screen parallel
         if _OPTIONS['with-screen-parallel'] then
            defines { "LSGL_ENABLE_SCREEN_PARALLEL" }
         end

         -- MPI
         if _OPTIONS['with-mpi'] then
            defines { "LSGL_ENABLE_MPI" }
         end

         if k_platform then 
            --buildoptions { "-Kfast", "-KPIC" }
            buildoptions { "-Xg -KPIC" }

            -- fj openmp
            if _OPTIONS['with-openmp'] then
               buildoptions { "-Kopenmp" }
            end
         else
            -- gcc openmp
            if _OPTIONS['with-openmp'] then
               buildoptions { "-fopenmp" }
               linkoptions { "-fopenmp" }
            end
         end

      configuration "Debug"
         defines { "DEBUG" } -- -DDEBUG
         if k_platform then 
            flags { "Symbols" }
         else
            flags { "Symbols", "EnableSSE2" }
         end

         if _OPTIONS['with-mpi'] then
            if _OPTIONS['with-screen-parallel'] then
               targetname "LSGLESd_mpi_screenparallel"
            else
               targetname "LSGLESd_mpi"
            end
         else
            targetname "LSGLESd"
         end

      configuration "Release"
         -- defines { "NDEBUG" } -- -NDEBUG

         if k_platform then 
            -- '-g' invoke -O0 optimization in FCCpx
            buildoptions { "-Kfast" }
         else
            flags { "Symbols", "Optimize", "EnableSSE2" }
         end

         if _OPTIONS['with-mpi'] then
            if _OPTIONS['with-screen-parallel'] then
               targetname "LSGLES_mpi_screen_parallel"
            else
               targetname "LSGLES_mpi"
            end
         else
            targetname "LSGLES"
         end

   if _OPTIONS['with-unittest'] then
      dofile "premake4-test.lua"
   end
