"""Independent mathematical oracle for the frozen production extrema expressions (offline).

The test compiles selected unchanged Fortran expressions from the actual source,
not a Python implementation of their bit packing. The oracle compares ordinary
decoded IEEE values and spatial tuple order. Huge formal dimensions need no
grid allocation in this scalar boundary projection.
"""
import argparse
import hashlib
import json
from pathlib import Path
import random
import struct
import subprocess as sp
import tempfile


def section(s, left, right):
    if s.count(left) != 1:
        raise ValueError('production anchor ambiguous: ' + left)
    after=s[s.index(left):]
    if right not in after:
        raise ValueError('missing production terminator: ' + right)
    return after[:after.index(right)]


def probe(src):
    code=(src/'baro_check.f').read_text()
    h=code.split('      RECURSIVE SUBROUTINE SMPI_TRY_HEALTH(',1)[1]
    prev=code.split('      RECURSIVE SUBROUTINE SMPI_VERDICT2(',1)[1].split(
        '      RECURSIVE SUBROUTINE SMPI_TRY_HEALTH(',1)[0]
    exact=lambda s:section(s,'        DO I=IMM1,1,-1','          IF(VMAXTMP.EQ.VAMAX)')
    assert exact(h)==exact(prev) and 'INTEGER, VOLATILE :: IPSPEED1' in h and 'IPSPEED1=I' in h
    assert 'IF(SMPI_OWNS(JVEL1)) THEN' in h and 'IF(VMAXTMP.EQ.VAMAX) THEN' in h
    rt=(src/'secom_mpi.f90').read_text()
    assert 'MPI_Allreduce(MPI_IN_PLACE,keys,3,MPI_INTEGER8,MPI_MAX' in rt
    def line(token):
        ls=[x for x in h.splitlines() if x.strip().startswith(token)]
        if len(ls)!=1:raise ValueError('ambiguous expression: '+token)
        return ls[0]+'\n'
    return ("""      PROGRAM PROBE
      USE ISO_FORTRAN_ENV, ONLY: INT64
      IMPLICIT NONE
      INTEGER IM,JM,I,J,BROW,BCELL,BSPEED,IERR,ISIZE
      INTEGER JW1,JVEL1,IAMAX,JAMAX,MW1,BW1,MV1,BV1,MZ1,BZ1
      INTEGER(INT64) KEYS1(3),BASE1,TOTAL1,SCORE1,TOP1,LOC1
      REAL VAMAX,ELMINX
""" + line('PARAMETER(BASE1=') + """
      DO
      READ(*,*,IOSTAT=IERR) IM,JM,I,J,BROW,BCELL,BSPEED
      IF(IERR.NE.0) EXIT
      KEYS1=0_INT64
      JW1=0
      JVEL1=0
      IAMAX=0
      JAMAX=0
      VAMAX=0.
      ELMINX=0.
""" + section(h,'      TOTAL1=INT(IM,INT64)','      DO J=SMPI_JLO,SMPI_JHI') + """
      BZ1=BROW
""" + section(h.split('      IF(JW1.NE.0) THEN',1)[0],
              '        MZ1=IBCLR(BZ1,31)','        IF(J.EQ.JM) CYCLE') + """
      DO ISIZE=1,1
        IF(J.EQ.JM) CYCLE
        BV1=BSPEED
""" + section(h,'        MV1=IBCLR(BV1,31)','      IF(JW1.NE.0) THEN') + """
      IF(JW1.NE.0) THEN
        MW1=IBCLR(BW1,31)
        BZ1=BCELL
        MZ1=IBCLR(BZ1,31)
        IF(BZ1.NE.BW1.AND..NOT.(MW1.EQ.0.AND.MZ1.EQ.0)) GOTO 71
""" + line('LOC1=INT(JW1') + line('KEYS1(2)=TOP1') + """
      ENDIF
      IF(KEYS1(1).NE.0_INT64) THEN
""" + section(h,'      MV1=INT(KEYS1(1)/BASE1)','      IF(SMPI_OWNS(JVEL1)) THEN') + """
      ENDIF
      IF(KEYS1(2).NE.0_INT64) THEN
""" + section(h,'      TOP1=KEYS1(2)/BASE1','      END\n') + """
      ENDIF
      GOTO 72
   71 KEYS1(3)=1
   72 WRITE(*,'(8(I0,1X))') KEYS1(1),KEYS1(2),KEYS1(3),
     & TRANSFER(VAMAX,0),JVEL1,TRANSFER(ELMINX,0),IAMAX,JAMAX
      ENDDO
      END
""")


def num(bits):
    return struct.unpack('>f',struct.pack('>I',bits))[0]


def signed(bits):
    return int.from_bytes(bits.to_bytes(4,'big'),'big',signed=True)


def run(src,fc):
    # Representatives of both ends and interior of every binary32 exponent bin.
    vals=[(sgn<<31)|(exp<<23)|mant for sgn in (0,1) for exp in range(1,255)
          for mant in (0,1,0x3fffff,0x7ffffe,0x7fffff)]
    vals += [0,1<<31,0x501502f8,0x501502f9,0x501502fa]
    seed=846313; rng=random.Random(seed)
    groups=[]; records=[]
    def add(im,jm,ix,j,w,c,v):
        records.append((im,jm,ix,j,w,c,v))
        return len(records)-1
    # One long group is a total-order comparison, not just a single minimum.
    for v in vals:
        add(32768,32768,32768,32767,v,v,0)
    for dims in ((372,479),(32768,32768),(2,1<<29),(1,(1<<30)-1),((1<<30),1)):
        im,jm=dims
        for w in (0,1<<31,0x00800000,0xff7fffff,0x7f7fffff,0x501502f8,0x501502f9,0x501502fa):
            for ix,j in ((1,1),(im,1),(1,jm),(im,jm),(im,max(jm-1,1))):
                add(im,jm,ix,j,w,(w^(1<<31)) if num(w)==0 else w,0x42c80000)
        for _ in range(256):
            # Distinct J: one speed row summary per owner. Ties are deliberate.
            js=rng.sample(range(1,jm+1),min(24,jm))
            same_w=rng.choice([0,0xbf800000,0x00800000,0xff7fffff,0x501502f8])
            same_v=rng.choice([0,0x42c80000,0x00800000])
            g=[]
            for j in js:
                w=rng.choice(vals) if rng.random()<.25 else same_w
                c=w^(1<<31) if num(w)==0 and rng.random()<.5 else w
                v=same_v^(1<<31) if same_v==0 and rng.random()<.5 else same_v
                g.append(add(im,jm,rng.randint(1,im),j,w,c,v))
            groups.append(g)
    # Targeted same-J/different-I water ties, with conflicting zero signs.
    for w in (0,1<<31,0xc4800000,0x501502f8,0xff7fffff):
        groups.append([add(32768,32768,i,1,w,w^(1<<31) if i==32768 and num(w)==0 else w,0)
                       for i in (1,2,16384,32768)])
    for im,jm in ((1,1<<30),(3,(1<<29)),(32769,32768)):
        add(im,jm,1,1,0,0,0)
    for v in (0x42c7ffff,0x42c80000,0x42c80001,0x7f7fffff,0xff7fffff,0x80000000):
        add(32768,32768,1,32767,0xbf800000,0xbf800000,v)
    # Exercise every normal speed exponent/mantissa witness, including veto >100.
    for v in vals:
        add(372,479,1,477,0xbfe00000,0xbfe00000,v)
    jg=[add(32768,32768,1,j,0,0,v) for j,v in [(1,0),(100,1<<31),(16000,0),(32767,1<<31)]]
    groups.append(jg)
    rows=''.join('%d %d %d %d %d %d %d\n'%(im,jm,i,j,signed(w),signed(c),signed(v))
                 for im,jm,i,j,w,c,v in records)
    code=probe(src)
    with tempfile.TemporaryDirectory() as tmp:
        tmp=Path(tmp); exe=tmp/'probe'
        sp.run([fc,'-x','f77','-ffixed-line-length-132','-std=legacy','-O0','-ftrapv','-fcheck=all',
                '-fno-fast-math','-ffpe-trap=invalid,zero,overflow','-','-o',str(exe)],
               input=code.encode(),check=True)
        out=sp.run([str(exe)],input=rows.encode(),stdout=sp.PIPE,check=True).stdout.splitlines()
    assert len(out)==len(records)
    ks=[];water_keys={};speed_keys={}
    for rid,(line,row) in enumerate(zip(out,records)):
        vk,wk,veto,vbits,vj,wbits,ix,jy=map(int,line.split())
        im,jm,i,j,w,c,v=row;speed=num(v);water=num(w)
        bad=(im*jm>1<<30 or jm>=1<<30 or (j!=jm and speed>100))
        assert bool(veto)==bad,(rid,row,'veto')
        ks.append((vk,wk))
        if bad:continue
        wvalid=water<1e10;vvalid=j!=jm and 0<=speed<=100
        assert 0<=vk<1<<63 and 0<=wk<1<<63
        assert bool(wk)==wvalid and bool(vk)==vvalid
        if wvalid:assert (ix,jy,wbits)==(i,j,signed(c))
        if vvalid:assert (vj,vbits)==(j,signed(v))
        # Same (J,I): isolate numerical ordering instead of comparing formulas.
        if rid<len(vals) and wvalid and water!=0:water_keys[water]=wk
        if im==372 and j==477 and vvalid and speed!=0:speed_keys[speed]=vk
    for d in (water_keys,speed_keys): assert len(d)>100
    assert all(water_keys[x]>water_keys[y] for x,y in zip(sorted(water_keys),sorted(water_keys)[1:]))
    assert all(speed_keys[x]<speed_keys[y] for x,y in zip(sorted(speed_keys),sorted(speed_keys)[1:]))
    assert max(wk for vk,wk in ks)==0x7fbfffffffffffff
    def wk_tuple(g):
        a=records[g];return (num(a[4]),a[3],a[2])
    def vk_tuple(g):
        a=records[g];return (num(a[6]),a[3])
    for group in groups:
        gw=[i for i in group if num(records[i][4])<1e10]
        if gw:
            want=min(gw,key=wk_tuple); got=max(group,key=lambda i:ks[i][1])
            assert wk_tuple(want)==wk_tuple(got) and records[want][5]==records[got][5]
        gv=[i for i in group if records[i][3]!=records[i][1]]
        if gv:
            want=max(gv,key=vk_tuple);got=max(group,key=lambda i:ks[i][0])
            assert vk_tuple(want)==vk_tuple(got) and records[want][6]==records[got][6]
    return {'pass':True,'seed':seed,'compiled_records':len(records),'numeric_order_groups':len(groups),
            'distinct_water_values_ordered':len(water_keys),'distinct_speed_values_ordered':len(speed_keys),
            'water_max_observed_legal_key':max(wk for vk,wk in ks),
            'signed_int64_max':(1<<63)-1,'sentinel_1e10_bits':'501502f9',
            'tested_total_cells_max':1<<30,'production_sha256':
            {n:hashlib.sha256((src/n).read_bytes()).hexdigest() for n in ('baro_check.f','secom_mpi.f90')},
            'compiled_projection_sha256':hashlib.sha256(code.encode()).hexdigest(),
            'method':'actual checked Fortran packing/decoding fragments vs independent float/J/I tuple oracle; not a solver run'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('source',type=Path)
    p.add_argument('--fc',default='gfortran');p.add_argument('--output',type=Path);a=p.parse_args()
    s=json.dumps(run(a.source,a.fc),indent=2)
    if a.output:
        with a.output.open('x') as f:f.write(s+'\n')
    print(s)
