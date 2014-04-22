for CMD in \
  "echo hello" \
  "echo fuck" \
  "echo love" \
  "echo dane" \
  "echo ahe" \
; do
  #echo -e "$ ${CMD}"
  eval "${CMD}"
done
