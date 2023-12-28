#!/usr/bin/perl
#May 9, 2012 by Terri Porter
#Perl script to conduct automated searches on BOLD's BIN database to retrieve species names
#usage perl search_BOLD.plx BIN.query
use strict;
use warnings;
# Create a user agent object
#use LWP::UserAgent;
use LWP::RobotUA;
use Web::Query 'wq';

#declare var
my $ua;
my $i=0;
my $binList;
my $j=0;
my $BIN;
my $req;
my $request;
my $res;
my $outfile;

#declare array
my @in;
my @binList;

$ua = LWP::RobotUA->new('SearchBINRobot/0.1','terriblue2002@yahoo.com');
#$ua = LWP::UserAgent->new;
#$ua->agent("MyApp/0.1 ");
$ua->delay(1); #max one page per minute

open (IN,"<",$ARGV[0]) || die "Error cannot open infile: $!\n";
@in = <IN>;
close IN;

while ($in[$i]) {
	$binList = $in[$i];
	chomp $binList;

	@binList = split(" ",$binList);

	while ($binList[$j]) {
		$BIN = $binList[$j];

		$req = "http://www.boldsystems.org/index.php/Public_BarcodeCluster?clusterguid=".$BIN;

		# Create a request
		$request = HTTP::Request->new(GET=> $req);

		# Pass request to the user agent and get a response back
		$res = $ua->request($request);

		$outfile = $BIN.".html";

		open (OUT,">>",$outfile) || die "Error cannot open outfile: $!\n";

		# Check the outcome of the response
		if ($res->is_success) {
			print OUT $res->content;
		}
		else {
			print OUT $res->status_line, "\n";
		}
		close OUT;
		$j++;
		$BIN=();
		$req=();
		$request=();
		$res=();
		$outfile=();
	}
	$j=0;
	$i++;
	$binList=();
	@binList=();
}
$i=0;

