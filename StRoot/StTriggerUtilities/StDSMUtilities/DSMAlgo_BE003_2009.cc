//
// Pibero Djawotho <pibero@comp.tamu.edu>
// Texas A&M University Cyclotron Institute
// 7 Jan 2009
//

#include "DSM.hh"
// #include "sumTriggerPatchChannels.hh"
#include "DSMAlgo_BE003_2009.hh"

void DSMAlgo_BE003_2009::operator()(DSM& dsm)
{
  // INPUT:

  // 10 12-bit BEMC channels
  // (5-0) high-tower
  // (11-6) trigger-patch

  // JP1 (0-15)  - ch1/3/5/7/9 (odd  channels - to upper DSM)
  // JP6 (16-31) - ch0/2/4/6/8 (even channels - to lower DSM)

  // REGISTERS:

  // R0: BEMC-High-Tower-Th0 (6)
  // R1: BEMC-High-Tower-Th1 (6)
  // R2: BEMC-High-Tower-Th2 (6)
  // R3: BEMC-High-Tower-Th3 (6)

  // ACTION:

  // JP1 - even channels - to lower DSM

  // int highTowerBitsJP1 = 0;
  // int  lowEtaSumJP1 = 0;
  // int highEtaSumJP1 = 0;

  // // Args: dsm, chMin, chMax, step, targetPedestal, sum, highTowerBits

  // sumTriggerPatchChannels(dsm, 0, 4, 2, 1,  lowEtaSumJP1, highTowerBitsJP1);
  // sumTriggerPatchChannels(dsm, 6, 8, 2, 1, highEtaSumJP1, highTowerBitsJP1);

  // // JP6 - odd channels - to upper DSM

  // int highTowerBitsJP6 = 0;
  // int  lowEtaSumJP6 = 0;
  // int highEtaSumJP6 = 0;

  // // Args: dsm, chMin, chMax, step, targetPedestal, sum, highTowerBits

  // sumTriggerPatchChannels(dsm, 1, 5, 2, 1,  lowEtaSumJP6, highTowerBitsJP6);
  // sumTriggerPatchChannels(dsm, 7, 9, 2, 1, highEtaSumJP6, highTowerBitsJP6);

  unsigned int highTowerBits[10][6];
  unsigned int trigPatchBits[10];

  for(int ichn = 0; ichn < 10; ichn++){
    unsigned int ht = dsm.channels[ichn] & 0x3f;
    for(int ireg = 0; ireg < 6; ireg++){
      highTowerBits[ichn][ireg] = ht > dsm.registers[ireg];
    }
    unsigned int tp = dsm.channels[ichn] >> 6 & 0x3f;
    trigPatchBits[ichn] = tp > dsm.registers[5];
  }

  unsigned int evenhtBits[6];
  unsigned int oddhtBits[6];

  for(int ireg = 0; ireg < 6; ireg++){
    evenhtBits[ireg] = 0;
    oddhtBits[ireg] = 0;

    for(int iichn = 0; iichn < 5; iichn++){
      evenhtBits[ireg] |= highTowerBits[2*iichn][ireg];
      oddhtBits[ireg] |= highTowerBits[2*iichn+1][ireg];
    }
  }

  unsigned int eventpBits = 0;
  unsigned int oddtpBits = 0;

  for(int iichn = 0; iichn < 5; iichn++){
    eventpBits |= trigPatchBits[2*iichn];
    oddtpBits |= trigPatchBits[2*iichn+1];
  }

  unsigned int evenhttpBits = 0;
  unsigned int oddhttpBits = 0;
  for(int iichn = 0; iichn < 5; iichn++){
    evenhttpBits |= (highTowerBits[2*iichn][5] & trigPatchBits[2*iichn]);
    oddhttpBits |= (highTowerBits[2*iichn+1][5] & trigPatchBits[2*iichn+1]);
  }


  // OUTPUT (32):

  // JP1 (0-15) to upper DSM

  // (0-5) TP sum for low-eta group (6)
  // (6-11) TP sum for high-eta group (6)
  // (12-15) HT bits (4)

  // JP6 (16-31) to lower DSM

  // (16-21) TP sum for low-eta group (6)
  // (22-27) TP sum for high-eta group (6)
  // (28-31) HT bits (4)

  // JP1 (0-15)

  // int out = 0;

  // out |=  lowEtaSumJP1;
  // out |= highEtaSumJP1    <<  6;
  // out |= highTowerBitsJP1 << 12;

  // // JP6 (16-31)

  // out |=  lowEtaSumJP6    << 16;
  // out |= highEtaSumJP6    << 22;
  // out |= highTowerBitsJP6 << 28;

  // OUTPUT (16):

  // (0-8) Unused
  // (9) TP threshold bit odd channels
  // (10) HT.TP threshold bit odd channels
  // (11-15) HT threshold bits odd channels
  // (16-24) Unused
  // (25) TP threshold bit even channels
  // (26) HT.TP threshold bit even channels
  // (27-31) HT threshold bit even channels
  int out = 0;

  out |=  oddtpBits << 8;
  out |=   oddhttpBits   <<  9;
  out |=   oddhtBits[0] << 10;
  out |=   oddhtBits[1] << 11;
  out |=   oddhtBits[2] << 12;
  out |=   oddhtBits[3] << 13;
  out |=   oddhtBits[4] << 14;
  out |=   oddhtBits[5] << 15;

  out |=  eventpBits << 24;
  out |=   evenhttpBits   <<  25;
  out |=   evenhtBits[0] << 26;
  out |=   evenhtBits[1] << 27;
  out |=   evenhtBits[2] << 28;
  out |=   evenhtBits[3] << 29;
  out |=   evenhtBits[4] << 30;
  out |=   evenhtBits[5] << 31;

  dsm.output = out;
}
