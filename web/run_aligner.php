<?php
	$sequences=$_POST['text'];
	`echo "$sequences" > input.fasta`;
	$gap_open=$_POST['gap_open'];
	$gap_extension=$_POST['gap_extension'];
	$algo=$_POST['algo'];
	//header('Location: loading.html', true);
	echo "<img src='loading.gif'>";
	`./multy -i input.fasta -g $gap_open -e $gap_extension $algo > out.txt`;
	header('Location: results.html', true);
?>