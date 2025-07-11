#!/bin/bash
declare -a nufilled=()

declare -a nuplaylistall=(
[0]=1A
[1]=1B
[2]=1C
[3]=1D
[4]=1E
[5]=1F
[6]=1G
[7]=1L
[8]=1M
[9]=1N
[10]=1O
[11]=1P
)

declare -a anufilled=()

declare -a anuplaylistall=(
[0]=6A
[1]=6B
[2]=6C
[3]=6D
[4]=6E
[5]=6F
[6]=6G
[7]=6H
[8]=6I
[9]=1J
)

declare -a materials=(
[0]=Lead
[1]=Iron
[2]=Carbon
)

declare -a runtypes=(
[0]=MC
[1]=Data
)


for dir in ${1}/*; do 
    if [ -d "$dir" ]; then 
        cd $dir
        hadd -f runEventLoopTargetsMCLead.root runEventLoopTargetsMC2082.root runEventLoopTargetsMC3082.root runEventLoopTargetsMC4082.root runEventLoopTargetsMC5082.root
        hadd -f runEventLoopTargetsMCIron.root runEventLoopTargetsMC2026.root runEventLoopTargetsMC3026.root  runEventLoopTargetsMC5026.root
        hadd -f runEventLoopTargetsMCCarbon.root runEventLoopTargetsMC3006.root
        hadd -f runEventLoopTargetsDataLead.root runEventLoopTargetsData2082.root runEventLoopTargetsData3082.root runEventLoopTargetsData4082.root runEventLoopTargetsData5082.root
        hadd -f runEventLoopTargetsDataIron.root runEventLoopTargetsData2026.root runEventLoopTargetsData3026.root  runEventLoopTargetsData5026.root
        hadd -f runEventLoopTargetsDataCarbon.root runEventLoopTargetsData3006.root
    fi 
done

combinednudir=${1}/combinedNu/
if [ ! -d ${combinednudir} ]; then
  mkdir -p ${combinednudir} 
fi

combinedanudir=${1}/combinedANu/
if [ ! -d ${combinedanudir} ]; then
  mkdir -p ${combinedanudir} 
fi


#Neutrino mode
#==========================
for rt in ${runtypes[@]}; do   
    for mat in ${materials[@]}; do   
        base="hadd -f ${combinednudir}/runEventLoopTargets${rt}${mat}.root"
        for dir in ${nuplaylistall[@]}; do   
            #echo -e $dir
            base="${base} ${1}/${dir}/runEventLoopTargets${rt}${mat}.root"
        done
        $base
    done
done
#Adding Water
for rt in ${runtypes[@]}; do   
    for mat in ${materials[@]}; do   
        base="hadd -f ${combinednudir}/runEventLoopTargets${rt}${mat}.root"
        for dir in ${nufilled[@]}; do   
            #echo -e $dir
            base="${base} ${1}/${dir}/runEventLoopTargets${rt}${mat}.root"
        done
        $base
    done
done
#anti-neutrino mode
#==========================
for rt in ${runtypes[@]}; do   
    for mat in ${materials[@]}; do   
        base="hadd -f ${combinednudir}/runEventLoopTargets${rt}${mat}.root"
        for dir in ${nuplaylistall[@]}; do   
            #echo -e $dir
            base="${base} ${1}/${dir}/runEventLoopTargets${rt}${mat}.root"
        done
        $base
    done
done
#Adding Water
for rt in ${runtypes[@]}; do   
    for mat in ${materials[@]}; do   
        base="hadd -f ${combinednudir}/runEventLoopTargets${rt}${mat}.root"
        for dir in ${anufilled[@]}; do   
            #echo -e $dir
            base="${base} ${1}/${dir}/runEventLoopTargets${rt}${mat}.root"
        done
        $base
    done
done