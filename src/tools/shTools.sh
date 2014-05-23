
function plural()
{
#return a word with or without 's' (plural) if the quantity is greater than 1
#call as plural(quantity, word)
#quantity (#1) is the quantity
#word (#2) is the word to put to plural
#
  if [ "$1" -gt 1 ]; then
    echo "$1 $2"s
  else
    echo "$1 $2"
  fi
}

function formatDateDiff(){
#format a number giving the difference between 2 dates, the result is given in terms of days, hours, minutes and second
#call as formatDateDiff number
#where number is the difference in seconds between the two dates

  dhms=$1
  d_day=$(( $dhms/(60*60*24) ))
  hms=$(( $dhms%(60*60*24) ))
  d_hour=$(( $hms/(60*60) ))
  ms=$(( $hms%(60*60) ))
  d_min=$(( $ms/60 ))
  d_sec=$(( $ms%60 ))
  #processing days
  if [ "${d_day}" -gt 0 ]; then
    duration=$( plural ${d_day} day )
  fi

  #processing hours
  if [ ! -z "${duration}" ]; then
    temp="${duration}"
    duration="${temp} $( plural ${d_hour} hour )"
  else
    if [ "${d_hour}" -gt 0 ]; then
      duration="$( plural ${d_hour} hour )"
    fi
  fi

  #processing minutes
  if [ ! -z "${duration}" ]; then
    temp="${duration}"
    duration="${temp} $( plural ${d_min} minute )"
  else
    if [ "${d_min}" -gt 0 ]; then
      duration="$( plural ${d_min}  minute )"
    fi
  fi
  #processing seconds
  if [ ! -z "${duration}" ]; then
    temp="${duration}"
    duration="${temp} $( plural ${d_sec} second )"
  else
    duration="$( plural ${d_sec}  second )"
  fi
  echo "${duration}"
}
