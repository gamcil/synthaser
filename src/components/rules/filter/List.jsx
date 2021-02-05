import React, { useMemo } from 'react'

import FilterItem from './Item'

const FilterList = ({
  filters,
  domains,
  allDomains, 
  handleAdd,
  handleChange,
  handleRemove,
}) => {

  const typeOptions = useMemo(
    () => domains.map(d => ({ label: d, value: d })),
    [domains]
  )

  const filterItems = filters.map(filter => (
    <FilterItem
      key={filter.uuid}
      type={filter.type}
      domains={filter.domains}
      allDomains={allDomains}
      typeOptions={typeOptions}
      handleChange={handleChange(filter.uuid)}
      handleRemove={handleRemove(filter.uuid)}
    />
  ))

  return (
    <div>
      <div>
        <button type="button" onClick={handleAdd}>Add</button>
      </div>
      <ul>
        {filterItems}
      </ul>
    </div>
  )
}

export default React.memo(FilterList)
