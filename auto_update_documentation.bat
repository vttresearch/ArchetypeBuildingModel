TITLE Auto update documentation
CALL git merge main
CALL julia auto_update_documentation.jl
CALL git commit -a -m "Update online documentation"
CALL git push