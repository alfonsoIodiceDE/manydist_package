# FIFA 21 Player Data - Dutch League

The `fifa_nl` dataset contains information on players in the Dutch
League from the FIFA 21 video game. This dataset includes various
attributes of players, such as demographics, club details, skill
ratings, and physical characteristics.

## Usage

``` r
data("fifa_nl")
```

## Format

A data frame with observations on various attributes describing the
players.

- `player_positions`:

  Primary playing positions of the player.

- `nationality`:

  The country the player represents.

- `team_position`:

  Player's assigned position within their club.

- `club_name`:

  Name of the club the player is part of.

- `work_rate`:

  The player's work rate, describing defensive and attacking intensity.

- `weak_foot`:

  Skill rating for the player's non-dominant foot, ranging from 1 to 5.

- `skill_moves`:

  Skill moves rating, indicating technical skill and ability to perform
  complex moves, on a scale of 1 to 5.

- `international_reputation`:

  Player's reputation on an international scale, from 1=local to
  3=global star.

- `body_type`:

  Body type of the player ( Lean, Normal, Stocky.

- `preferred_foot`:

  Dominant foot of the player, either Left or Right.

- `age`:

  Age of the player in years.

- `height_cm`:

  Height of the player in centimeters.

- `weight_kg`:

  Weight of the player in kilograms.

- `overall`:

  Overall skill rating of the player out of 100.

- `potential`:

  Potential skill rating the player may achieve in the future.

- `value_eur`:

  Estimated market value of the player in Euros.

- `wage_eur`:

  Player's weekly wage in Euros.

- `release_clause_eur`:

  Release clause value in Euros, which other clubs must pay to buy out
  the player's contract.

- `pace`:

  Speed rating of the player out of 100.

- `shooting`:

  Shooting skill rating out of 100.

- `passing`:

  Passing skill rating out of 100.

- `dribbling`:

  Dribbling skill rating out of 100.

- `defending`:

  Defending skill rating out of 100.

- `physic`:

  Physicality rating out of 100.

## Details

This dataset provides a snapshot of player attributes and performance
indicators as represented in FIFA 21 for players in the Dutch League. It
can be used to analyze player characteristics, compare skills across
players, and explore potential relationships among variables such as
age, position, and various skill ratings.

## References

Stefano Leone. (2021). *FIFA 21 Complete Player Dataset*. Retrieved from
<https://www.kaggle.com/datasets/stefanoleone992/fifa-21-complete-player-dataset>.

## Examples

``` r
data(fifa_nl)
summary(fifa_nl)
#>  player_positions      nationality  team_position        club_name  
#>  CB     : 58      Netherlands:220   SUB    :186   ADO Den Haag: 27  
#>  ST     : 39      Germany    : 26   RES    : 42   FC Emmen    : 27  
#>  LB     : 26      Belgium    : 20   LCB    : 18   RKC Waalwijk: 26  
#>  RB     : 26      Sweden     : 10   RCB    : 18   Willem II   : 25  
#>  CDM, CM: 18      Morocco    :  9   LB     : 17   FC Groningen: 24  
#>  CM, CDM: 17      Norway     :  9   RB     : 17   FC Utrecht  : 24  
#>  (Other):224      (Other)    :114   (Other):110   (Other)     :255  
#>          work_rate   weak_foot skill_moves international_reputation
#>  Medium/Medium:202   2: 67     2:186       1:380                   
#>  High/Medium  : 77   3:264     3:187       2: 21                   
#>  Medium/High  : 46   4: 69     4: 31       3:  7                   
#>  Medium/Low   : 29   5:  8     5:  4                               
#>  High/High    : 23                                                 
#>  High/Low     : 20                                                 
#>  (Other)      : 11                                                 
#>   body_type   preferred_foot      age          height_cm       weight_kg    
#>  Lean  :115   Left :123      Min.   :17.00   Min.   :166.0   Min.   :55.00  
#>  Normal:272   Right:285      1st Qu.:20.00   1st Qu.:176.0   1st Qu.:69.00  
#>  Stocky: 21                  Median :22.00   Median :180.0   Median :74.00  
#>                              Mean   :23.27   Mean   :180.7   Mean   :73.59  
#>                              3rd Qu.:26.00   3rd Qu.:185.0   3rd Qu.:78.00  
#>                              Max.   :36.00   Max.   :201.0   Max.   :94.00  
#>                                                                             
#>     overall        potential       value_eur           wage_eur    
#>  Min.   :53.00   Min.   :62.00   Min.   :  120000   Min.   :  500  
#>  1st Qu.:63.00   1st Qu.:69.00   1st Qu.:  475000   1st Qu.: 2000  
#>  Median :66.00   Median :73.00   Median :  800000   Median : 3000  
#>  Mean   :66.53   Mean   :73.37   Mean   : 2343652   Mean   : 4768  
#>  3rd Qu.:70.00   3rd Qu.:76.25   3rd Qu.: 1900000   3rd Qu.: 6000  
#>  Max.   :84.00   Max.   :88.00   Max.   :27500000   Max.   :30000  
#>                                                                    
#>  release_clause_eur      pace          shooting        passing     
#>  Min.   :       0   Min.   :31.00   Min.   :24.00   Min.   :29.00  
#>  1st Qu.:  683000   1st Qu.:63.00   1st Qu.:44.00   1st Qu.:53.00  
#>  Median : 1100000   Median :69.00   Median :56.00   Median :59.00  
#>  Mean   : 3547203   Mean   :68.88   Mean   :53.83   Mean   :58.78  
#>  3rd Qu.: 2725000   3rd Qu.:76.00   3rd Qu.:64.00   3rd Qu.:65.00  
#>  Max.   :38500000   Max.   :93.00   Max.   :82.00   Max.   :84.00  
#>                                                                    
#>    dribbling       defending         physic     
#>  Min.   :35.00   Min.   :19.00   Min.   :36.00  
#>  1st Qu.:61.00   1st Qu.:36.00   1st Qu.:58.00  
#>  Median :65.00   Median :56.00   Median :66.00  
#>  Mean   :64.58   Mean   :51.17   Mean   :64.13  
#>  3rd Qu.:70.00   3rd Qu.:63.25   3rd Qu.:72.00  
#>  Max.   :86.00   Max.   :83.00   Max.   :89.00  
#>                                                 
```
