//
// Pibero Djawotho <pibero@comp.tamu.edu>
// Texas A&M University Cyclotron Institute
// 7 Jan 2009
//

#include <algorithm>
#include "bits.hh"
#include "DSM.hh"
#include "DSMAlgo_EM201_2009.hh"

int DSMAlgo_EM201_2009::ajpBarrel(const DSM& dsm, int offset) const
{
  int jpBits[6];

  // BC101-106

  for (int ch = 0; ch < 6; ++ch)
    {
      jpBits[ch] = dsm.channels[ch] >> offset & 0x3;
//      printf("The channel %d jp bit is %d\n", ch, jpBits[ch]); //Test by Z. Chang
    }
  const int R3 = dsm.registers[3];
//  printf("R3 is %d out of %d\n", R3,dsm.registers[3]);

  return (((jpBits[0] > R3) && (jpBits[1] > R3)) ||
	  ((jpBits[1] > R3) && (jpBits[2] > R3)) ||
	  ((jpBits[2] > R3) && (jpBits[3] > R3)) ||
	  ((jpBits[3] > R3) && (jpBits[4] > R3)) ||
	  ((jpBits[4] > R3) && (jpBits[5] > R3)) ||
	  ((jpBits[5] > R3) && (jpBits[0] > R3)));
}

int DSMAlgo_EM201_2009::ajpEndcap(const DSM& dsm) const
{
  int jpBits[6];
  const int R3 = dsm.registers[3];

  // EE101

  jpBits[0] = dsm.channels[6]      & 0x3; // JPA (4 o'clock)
  jpBits[1] = dsm.channels[6] >> 2 & 0x3; // JPB (6 o'clock)
  jpBits[2] = dsm.channels[6] >> 4 & 0x3; // JPC (8 o'clock)

  // EE102

  jpBits[3] = dsm.channels[7]      & 0x3; // JPA (10 o'clock)
  jpBits[4] = dsm.channels[7] >> 2 & 0x3; // JPB (12 o'clock)
  jpBits[5] = dsm.channels[7] >> 4 & 0x3; // JPC (2  o'clock)

  return (((jpBits[0] > R3) && (jpBits[1] > R3)) ||
	  ((jpBits[1] > R3) && (jpBits[2] > R3)) ||
	  ((jpBits[2] > R3) && (jpBits[3] > R3)) ||
	  ((jpBits[3] > R3) && (jpBits[4] > R3)) ||
	  ((jpBits[4] > R3) && (jpBits[5] > R3)) ||
	  ((jpBits[5] > R3) && (jpBits[0] > R3)));
}

void DSMAlgo_EM201_2009::operator()(DSM& dsm)
{
  // INPUT:

  // EM201 - ch0 - BEMC BC101 - 10' - JP0 & JP6 (West & East)
  //         ch1 - BEMC BC102 - 12' - JP1 & JP7
  //         ch2 - BEMC BC103 -  2' - JP2 & JP8
  //         ch3 - BEMC BC104 -  4' - JP3 & JP9
  //         ch4 - BEMC BC105 -  6' - JP4 & JP10
  //         ch5 - BEMC BC106 -  8' - JP5 & JP11
  //         ch6 - EEMC EE101 -  4',  6' and 8' - JP3, JP4 & JP5
  //         ch7 - EEMC EE102 - 10', 12' and 2' - JP0, JP1 & JP2

  // From BC101-106 (16):

  // (0-1) JPX (east, -1 < eta < 0) threshold bits (2)
  // (2-3) JPY (middle, -0.6 < eta < 0.4) threshold bits (2)
  // (4-5) JPZ (west, 0 < eta < 1) threshold bits (2)
  // (6-11) JPpartial (0.4 < eta < 1) sum (6)
  // (12-15) HT bits (4)

  // From EE101 (16):

  // (0-1) JPA (4 o'clock) threshold bits (2)
  // (2-3) JPB (6 o'clock) threshold bits (2)
  // (4-5) JPC (8 o'clock) threshold bits (2)
  // (6-11) Selected partial jet patch sum (6)
  // (12-13) Partial jet patch ID (1=A, 2=B, 3=C) (2)
  // (14-15) HT bits (2)

  // From EE102 (16):

  // (0-1) JPA (10 o'clock) threshold bits (2)
  // (2-3) JPB (12 o'clock) threshold bits (2)
  // (4-5) JPC (2  o'clock) threshold bits (2)
  // (6-11) Selected partial jet patch sum (6)
  // (12-13) Partial jet patch ID (1=A, 2=B, 3=C) (2)
  // (14-15) HT bits (2)

  // REGISTERS:

  // R0: Hybrid jet patch threshold-0
  // R1: Hybrid jet patch threshold-1
  // R2: Hybrid jet patch threshold-2

  // ACTION:

  // Complete hybrid jet patches using partial jet patch ID from EEMC

  int jpSum1 = dsm.channels[6] >> 6 & 0x3f; // Partial sum from EE101
  int jpSum2 = dsm.channels[7] >> 6 & 0x3f; // Partial sum from EE102

  int jpId1 = dsm.channels[6] >> 12 & 0x3; // Partial jet patch ID from EE101
  int jpId2 = dsm.channels[7] >> 12 & 0x3; // Partial jet patch ID from EE102

  switch (jpId1) {
  case 1: jpSum1 += dsm.channels[3] >> 6 & 0x3f; break; // Add partial sum from BC104 (4')
  case 2: jpSum1 += dsm.channels[4] >> 6 & 0x3f; break; // Add partial sum from BC105 (6')
  case 3: jpSum1 += dsm.channels[5] >> 6 & 0x3f; break; // Add partial sum from BC106 (8')
  }

  switch (jpId2) {
  case 1: jpSum2 += dsm.channels[0] >> 6 & 0x3f; break; // Add partial sum from BC101 (10')
  case 2: jpSum2 += dsm.channels[1] >> 6 & 0x3f; break; // Add partial sum from BC102 (12')
  case 3: jpSum2 += dsm.channels[2] >> 6 & 0x3f; break; // Add partial sum from BC103 (2')
  }

  // Combine (OR) the HT bits from the six BEMC layer 1 DSM's

  int htBitsBarrel = 0;

  for (int ch = 0; ch < 6; ++ch)
    htBitsBarrel |= dsm.channels[ch] >> 12 & 0xf;

  // Combine (OR) the HT bits from the two EEMC layer 1 DSM's

  int htBitsEndcap = 0;

  for (int ch = 6; ch < 8; ++ch)
    htBitsEndcap |= dsm.channels[ch] >> 14 & 0x3;

  // Combine (OR) the JP bits for the BEMC and EEMC separately

  int jpBitsBarrel = 0;

  for (int ch = 0; ch < 6; ++ch) {
    int jpx = dsm.channels[ch]      & 0x3;
    int jpy = dsm.channels[ch] >> 2 & 0x3;
    int jpz = dsm.channels[ch] >> 4 & 0x3;

    if (jpx > jpBitsBarrel) jpBitsBarrel = jpx;
    if (jpy > jpBitsBarrel) jpBitsBarrel = jpy;
    if (jpz > jpBitsBarrel) jpBitsBarrel = jpz;
  }

  int bjp1 = jpBitsBarrel > 1;
  int bjp2 = jpBitsBarrel > 2;

  int jpBitsEndcap = 0;

  for (int ch = 6; ch < 8; ++ch) {
    int jpa = dsm.channels[ch]      & 0x3;
    int jpb = dsm.channels[ch] >> 2 & 0x3;
    int jpc = dsm.channels[ch] >> 4 & 0x3;

    if (jpa > jpBitsEndcap) jpBitsEndcap = jpa;
    if (jpb > jpBitsEndcap) jpBitsEndcap = jpb;
    if (jpc > jpBitsEndcap) jpBitsEndcap = jpc;
  }

  int ejp1 = jpBitsEndcap > 1;
  int ejp2 = jpBitsEndcap > 2;

  // Compare the two completed hybrid jet patches to three thresholds
  // and combine (OR) the results with the BEMC-only and EEMC-only bits

  int jpBits = 0;

  for (int reg = 0; reg < 3; ++reg)
    if (jpSum1 > dsm.registers[reg] || jpSum2 > dsm.registers[reg]) ++jpBits;

  if (jpBitsBarrel > jpBits) jpBits = jpBitsBarrel;
  if (jpBitsEndcap > jpBits) jpBits = jpBitsEndcap;

  int jp0 = jpBits > 0;
  int jp1 = jpBits > 1;
  int jp2 = jpBits > 2;

  // Adjacent jet patch logic

  int ajpx = ajpBarrel(dsm, 0);
  int ajpy = ajpBarrel(dsm, 2);
  int ajpz = ajpBarrel(dsm, 4);
  int bajp = ajpx || ajpy || ajpz;
  int eajp = ajpEndcap(dsm);
  int  ajp = bajp || eajp;

  // OUTPUT (16):

  // (0:3) Barrel HT bits (4)
  // (4:5) Endcap HT bits (2)
  // (6) JP1, unified over the BEMC+EEMC (1)
  // (7) JP2, unified over the BEMC+EEMC (1)
  // (8) BJP1 for the 18 BEMC-only patches (1)
  // (9) BJP2 for the 18 BEMC-only patches (1)
  // (10) EJP1 for the 6 EEMC-only patches (1)
  // (11) EJP2 for the 6 EEMC-only patches (1)
  // (12) AJP for BEMC and EEMC but NOT the boundary (1)
  // (13) BAJP for the BEMC-only patches (1)
  // (14) EAJP for the EEMC-only patches (1)
  // (15) JP0, unified over the BEMC+EEMC (1)

  int out = 0;

  out |= htBitsBarrel;
  out |= htBitsEndcap << 4;
  out |= jp1  << 6;
  out |= jp2  << 7;
  out |= bjp1 << 8;
  out |= bjp2 << 9;
  out |= ejp1 << 10;
  out |= ejp2 << 11;
  out |= ajp  << 12;
  out |= bajp << 13;
  out |= eajp << 14;
  out |= jp0  << 15;

  // // INPUT:

  // // 6 channels from BEMC

  // // BEMC BC101 10 o'clock
  // // BEMC BC102 12 o'clock
  // // BEMC BC103 2 o'clock
  // // BEMC BC104 4 o'clock
  // // BEMC BC105 6 o'clock
  // // BEMC BC106 8 o'clock

  // // bits 0-7 unused
  // // bit 8 DAQ10k test bit
  // // bit 9 TP threshold bit
  // // bit 10 HT.TP threshold bit
  // // bits 11-15 HT threshold bits

  // // 2 channels from EEMC

  // // EEMC EE101 4, 6 and 8 o'clock
  // // EEMC EE102 10, 12 and 2 o'clock

  // // bits 0-13 unused
  // // bits 14-15 HT threshold bits

  // // REGISTERS:

  // // R0: DAQ10k-Sector-Count(3)
  // // R1: EMC-UPC-Topo-Swith (3)
  // // ACTION:
  // const int R0 = dsm.registers[0];
  // const int R1 = dsm.registers[1];

  // int bemcHT = 0; 
  // int bemcTP = 0;
  // int bemcHTTP = 0;

  // // Combine (OR) the HT bits from the six BEMC layer 1 DSM's
  // int bemcTPBit[6];
  // int bemcHTTPBit[6];
  // int bemcHTUPCBit[6];

  // int counterDAQ10K = 0;

  // for (int ichn = 0; ichn < 6; ++ichn) {
  //   bemcHT |= dsm.channels[ichn] >> 10 & 0x3f;

  //   bemcTPBit[ichn] = dsm.channels[ichn] >> 8 & 0x1;
  //   bemcTP |= bemcTPBit[ichn];

  //   bemcHTTPBit[ichn] = dsm.channels[ichn] >> 9 & 0x1;
  //   bemcHTTP |= bemcHTTPBit[ichn];

  //   bemcHTUPCBit[ichn] = dsm.channels[ichn] >> 15 & 0x1;

  //   counterDAQ10K += dsm.channels[ichn] >> 7 & 0x1;
  // }
  // int bemcDAQ10K = counterDAQ10K >= R0;

  // // Combine (OR) the HT bits from the two EEMC layer 1 DSM's

  // int eemcHT = 0; 

  // for (int ichn = 6; ichn < 8; ++ichn) {
  //   eemcHT |= dsm.channels[ichn] >> 14 & 0x3;
  // }

  // int tpB2B = 0;
  // int tpNONADJ = 0;

  // int httpB2B = 0;
  // int httpNONADJ = 0;

  // int htUPCB2B = 0;
  // int htUPCNONADJ = 0;

  // for(int ichn = 0; ichn < 3; ichn++){
  //   int jchn = (ichn+3)%12;
  //   tpB2B |= bemcTPBit[ichn] && bemcTPBit[jchn];
  //   httpB2B |= bemcHTTPBit[ichn] && bemcHTTPBit[jchn];
  //   htUPCB2B |= bemcHTUPCBit[ichn] && bemcHTUPCBit[jchn];
  // }

  // tpNONADJ = (bemcTPBit[0] && (bemcTPBit[2] || bemcTPBit[3] || bemcTPBit[4]))
  //   || (bemcTPBit[1] && (bemcTPBit[3] || bemcTPBit[4] || bemcTPBit[5]))
  //   || (bemcTPBit[2] && (bemcTPBit[4] || bemcTPBit[5]))
  //   || (bemcTPBit[3] && bemcTPBit[5]);

  // httpNONADJ = (bemcHTTPBit[0] && (bemcHTTPBit[2] || bemcHTTPBit[3] || bemcHTTPBit[4]))
  //   || (bemcHTTPBit[1] && (bemcHTTPBit[3] || bemcHTTPBit[4] || bemcHTTPBit[5]))
  //   || (bemcHTTPBit[2] && (bemcHTTPBit[4] || bemcHTTPBit[5]))
  //   || (bemcHTTPBit[3] && bemcHTTPBit[5]);

  // htUPCNONADJ = (bemcHTUPCBit[0] && (bemcHTUPCBit[2] || bemcHTUPCBit[3] || bemcHTUPCBit[4]))
  //   || (bemcHTUPCBit[1] && (bemcHTUPCBit[3] || bemcHTUPCBit[4] || bemcHTUPCBit[5]))
  //   || (bemcHTUPCBit[2] && (bemcHTUPCBit[4] || bemcHTUPCBit[5]))
  //   || (bemcHTUPCBit[3] && bemcHTUPCBit[5]);

  // int tpTopo = 0;
  // if(btest(R1, 0)) tpTopo = tpB2B;
  // else tpTopo = tpNONADJ;

  // int httpTopo = 0;
  // if(btest(R1, 1)) httpTopo = httpB2B;
  // else httpTopo = httpNONADJ;

  // int htUPCTopo = 0;
  // if(btest(R1, 2)) htUPCTopo = htUPCB2B;
  // else htUPCTopo = htUPCNONADJ;
  // // OUTPUT (16):

  // // (0:4) Barrel HT bits (5)
  // // (5) Barrel HT UPC
  // // (6) Barrel TP bit
  // // (7) Barrel HT.TP bit
  // // (8) Barrel TP-based topo bit
  // // (9) Barrel HTTP-based topo bit
  // // (10) Barrel UPC-based topo bit
  // // (11) Unused
  // // (12) Unused
  // // (13:14) Endcap HT bits (2)
  // // (15) DAQ10k

  // int out = 0;

  // out |= (bemcHT & 0x3f);
  // out |= bemcTP << 6;
  // out |= bemcHTTP  << 7;
  // out |= tpTopo  << 8;
  // out |= httpTopo << 9;
  // out |= htUPCTopo << 10;
  // out |= eemcHT << 13;
  // out |= bemcDAQ10K << 15;


  dsm.output = out;

  // // INFO

  dsm.info[0] = jpSum1;
  dsm.info[1] = jpSum2;
}