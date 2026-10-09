MODULE PARAMETERS
Implicit None

! numero massimo di atomi
Integer, Parameter :: natmax=2000
! maximum numbre of neighbours
Integer, Parameter :: nvmax0=natmax
! maximum number of elements
Integer, Parameter :: nelemmax=10
! maximum number of files to be minimized
Integer, Parameter :: nmaxfiles=10000
! maximum number of characters in the name of an input file
Integer, Parameter :: ncharmax=40
! parametro di tolleranza sulla determinazione di un minimo locale
Real(8), Parameter :: ftoll=1.e-7_8
Real(8), PARAMETER :: pi=3.141592653589793
Real(8), PARAMETER :: cbol=8.62d-05
Real(8), PARAMETER :: evsujoule=1.60219d-19
Real(8), PARAMETER :: uasukg=1./6.023d26
Real(8), PARAMETER :: angsum=1.d-10
!names of the order parameters
END MODULE PARAMETERS
