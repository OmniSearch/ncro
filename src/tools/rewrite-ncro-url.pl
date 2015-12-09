my %registry;
my %prefixes;
my $counter=0;
my $seen;
  
sub seenOrRegister 
{ $seen = $registry{@_[0]};
  if ($seen) { return $seen }
  $id = @_[1];
  $new = sprintf("http://purl.obolibrary.org/obo/NCRO_%07d",$counter++);
  $registry{@_[0]}=$new;
  $prefix = @_[1];
  if ($prefix =~ /([A-Z]+)\d{3-7}/)
    { $prefixes{$1}++;}
  return $new
}
      
open INPUT, "</Users/lori/repos/ncro/src/ontology/ncro.owl";
while ($line=<INPUT>)
  { print $line;
    last if ($line =~ m|\s*</owl:Ontology>|)
  }
    
while ($line=<INPUT>)
  { $line=~ s/^(.*?)ncro_(\d+.*)/$1NCRO_$2/;
    $line=~s|"(http://purl.obolibrary.org/obo/ncro-obo-WorkingVersion-06302015#(.*?))"|seenOrRegister($1,$2)|e;
    print $line;
#    $DB::single=1
#  print $line;
}

$DB::single=1;
$DB::single=1;
      
