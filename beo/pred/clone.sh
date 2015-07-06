# change sim1 to simX
X="3"

for i in `seq 1 20`;
do
	echo $i

	FOUT="sim${X}_${i}.R"
	cat sim1_$i.R | sed "s/which_exp <- 1/which_exp <- ${X}/" > $FOUT

	rm sim1_$i.R
done

FOUT="submit${X}.sh"
cat submit1.sh | sed "s/sim1/sim${X}/g" > $FOUT
chmod +x $FOUT
rm submit1.sh
