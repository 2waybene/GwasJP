## This script will create all of your analysis directories in a four loop

declare -a dirs=(
		"icaps/icaps.thiazide.all_races"
		"icaps/icaps.thiazide.white"
		"icaps/icaps.thiazide.black"
		"icaps/icaps.thiazide.hispanic"

		"icaps/icaps.ccb.all_races"
		"icaps/icaps.ccb.white"
		"icaps/icaps.ccb.black"
		"icaps/icaps.ccb.hispanic"

		"icaps/icaps.bb.all_races"
		"icaps/icaps.bb.white"
		"icaps/icaps.bb.black"
		"icaps/icaps.bb.hispanic"

		"icaps/icaps.aceia.arb.all_races"
		"icaps/icaps.aceia.arb.white"
		"icaps/icaps.aceia.arb.black"
		"icaps/icaps.aceia.arb.hispanic"
		)

echo "mkdir ../icaps"
for i in "${dirs[@]}"
do
   mkdir "../$i"
done
