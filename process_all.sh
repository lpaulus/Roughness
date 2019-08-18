while read name; do
    echo "Processing $name with arguments $@"
    bash process.sh "$name" "$@"
done < all.txt
