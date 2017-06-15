# Install script for directory: /afs/ihep.ac.cn/users/e/evan/myWorkSpace/Time

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/afs/ihep.ac.cn/users/e/evan/myWorkSpace/Time")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "RelWithDebInfo")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

# Install shared libraries without execute permission?
IF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  SET(CMAKE_INSTALL_SO_NO_EXE "0")
ENDIF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FOREACH(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libRanger.so.3.2.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libRanger.so.3.2"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libRanger.so"
      )
    IF(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      FILE(RPATH_CHECK
           FILE "${file}"
           RPATH "/afs/ihep.ac.cn/users/e/evan/myWorkSpace/Time/lib:/afs/ihep.ac.cn/users/m/manqi/ArborV3Deve/SL6/Marlin/v01-05/lib:/afs/ihep.ac.cn/soft/common/gcc/v01-17-05/lcio/v02-04-03/lib:/afs/ihep.ac.cn/soft/common/gcc/v01-17-05/gear/v01-04/lib:/afs/ihep.ac.cn/soft/common/gcc/v01-17-05/CLHEP/2.1.3.1/lib:/afs/ihep.ac.cn/soft/common/gcc/v01-17-05/ilcutil/v01-01/lib:/afs/ihep.ac.cn/soft/common/gcc/v01-17-05/root/5.34.07/lib")
    ENDIF()
  ENDFOREACH()
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE FILES
    "/afs/ihep.ac.cn/users/e/evan/myWorkSpace/Time/build/lib/libRanger.so.3.2.1"
    "/afs/ihep.ac.cn/users/e/evan/myWorkSpace/Time/build/lib/libRanger.so.3.2"
    "/afs/ihep.ac.cn/users/e/evan/myWorkSpace/Time/build/lib/libRanger.so"
    )
  FOREACH(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libRanger.so.3.2.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libRanger.so.3.2"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libRanger.so"
      )
    IF(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      FILE(RPATH_CHANGE
           FILE "${file}"
           OLD_RPATH "/afs/ihep.ac.cn/users/m/manqi/ArborV3Deve/SL6/Marlin/v01-05/lib:/afs/ihep.ac.cn/soft/common/gcc/v01-17-05/lcio/v02-04-03/lib:/afs/ihep.ac.cn/soft/common/gcc/v01-17-05/gear/v01-04/lib:/afs/ihep.ac.cn/soft/common/gcc/v01-17-05/CLHEP/2.1.3.1/lib:/afs/ihep.ac.cn/soft/common/gcc/v01-17-05/ilcutil/v01-01/lib:/afs/ihep.ac.cn/soft/common/gcc/v01-17-05/root/5.34.07/lib::::::::::::::::::::::::::::::::::::::::::::::::::"
           NEW_RPATH "/afs/ihep.ac.cn/users/e/evan/myWorkSpace/Time/lib:/afs/ihep.ac.cn/users/m/manqi/ArborV3Deve/SL6/Marlin/v01-05/lib:/afs/ihep.ac.cn/soft/common/gcc/v01-17-05/lcio/v02-04-03/lib:/afs/ihep.ac.cn/soft/common/gcc/v01-17-05/gear/v01-04/lib:/afs/ihep.ac.cn/soft/common/gcc/v01-17-05/CLHEP/2.1.3.1/lib:/afs/ihep.ac.cn/soft/common/gcc/v01-17-05/ilcutil/v01-01/lib:/afs/ihep.ac.cn/soft/common/gcc/v01-17-05/root/5.34.07/lib")
      IF(CMAKE_INSTALL_DO_STRIP)
        EXECUTE_PROCESS(COMMAND "/afs/ihep.ac.cn/soft/common/gcc/v01-17-05/mysql/usr/bin/strip" "${file}")
      ENDIF(CMAKE_INSTALL_DO_STRIP)
    ENDIF()
  ENDFOREACH()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(CMAKE_INSTALL_COMPONENT)
  SET(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
ELSE(CMAKE_INSTALL_COMPONENT)
  SET(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
ENDIF(CMAKE_INSTALL_COMPONENT)

FILE(WRITE "/afs/ihep.ac.cn/users/e/evan/myWorkSpace/Time/build/${CMAKE_INSTALL_MANIFEST}" "")
FOREACH(file ${CMAKE_INSTALL_MANIFEST_FILES})
  FILE(APPEND "/afs/ihep.ac.cn/users/e/evan/myWorkSpace/Time/build/${CMAKE_INSTALL_MANIFEST}" "${file}\n")
ENDFOREACH(file)
