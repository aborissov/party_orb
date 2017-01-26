SDF_PATH='..IDL/SDF/IDL'
LARE_PATH='./IDL'
!path = SDF_PATH + PATH_SEP(/SEARCH_PATH) + !path
!path = LARE_PATH + PATH_SEP(/SEARCH_PATH) + !path

@./IDL/SDF/IDL/Start.pro
.r getenergy
.r getprobe
