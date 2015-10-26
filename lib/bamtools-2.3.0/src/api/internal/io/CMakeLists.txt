# ==========================
# BamTools CMakeLists.txt
# (c) 2011 Derek Barnett
#
# src/api/internal/io
# ==========================

set( InternalIODir "${InternalDir}/io" )

#--------------------------
# platform-independent IO
#--------------------------
set( CommonIOSources
        ${InternalIODir}/BamDeviceFactory_p.cpp
        ${InternalIODir}/BamFile_p.cpp
        ${InternalIODir}/BamFtp_p.cpp
        ${InternalIODir}/BamHttp_p.cpp
        ${InternalIODir}/BamPipe_p.cpp
        ${InternalIODir}/BgzfStream_p.cpp
        ${InternalIODir}/ByteArray_p.cpp
        ${InternalIODir}/HostAddress_p.cpp
        ${InternalIODir}/HostInfo_p.cpp
        ${InternalIODir}/HttpHeader_p.cpp
        ${InternalIODir}/ILocalIODevice_p.cpp
        ${InternalIODir}/RollingBuffer_p.cpp
        ${InternalIODir}/TcpSocket_p.cpp
        ${InternalIODir}/TcpSocketEngine_p.cpp
)

#------------------------
# platform-dependent IO
#------------------------
if( WIN32 )
    set( PlatformIOSources ${InternalIODir}/TcpSocketEngine_win_p.cpp )
else()
    set( PlatformIOSources ${InternalIODir}/TcpSocketEngine_unix_p.cpp )
endif()

#---------------------------
# make build-specific list
#---------------------------
set( InternalIOSources
        ${CommonIOSources}
        ${PlatformIOSources} 

        PARENT_SCOPE # <-- leave this last
)

