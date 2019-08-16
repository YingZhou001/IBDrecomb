../src/IBDrecomb -refinedibd test.ibd.gz -fbin 10000 > test.10000.map

../src/IBDrecomb -ibd test.2.ibd.gz -fbin 10000 > test.10000.map

../src/IBDrecomb -refinedibd test.ibd.gz -fbin 10000 -rate 1.21 > test.10000.normalized.map
